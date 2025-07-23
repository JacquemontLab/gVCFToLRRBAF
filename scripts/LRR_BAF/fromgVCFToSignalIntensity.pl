#!/usr/bin/env perl
use warnings;     # Display helpful warnings for potential issues
use strict;       # Enforce variable declarations to avoid bugs
use Getopt::Long; # Module to parse command-line options (like --sample_id)



# ------------------------------------------
# Script purpose:
# This script processes a compressed gVCF file (.gvcf.gz) to extract biallelic SNPs 
# and compute two types of variant-level metrics used in CNV and allele-specific analyses:
#
# 1. Coverage-based SNP metrics:
#    - Output file: <prefix>.snp_metrics.tsv
#    - Content: Per-SNP coverage (DP), B Allele Count (BAC, i.e., count of ALT alleles), 
#      and B Allele Frequency (BAF = BAC / DP)
#
# 2. Log R Ratio (LRR) and BAF for CNV detection tools:
#    - Output file: <prefix>.baf_lrr.tsv
#    - Content: Formatted for compatibility with CNV calling tools such as PennCNV and QuantiSNP.
#      Includes genomic coordinates, BAF, and LRR (log ratio of observed vs mean coverage).
#
# Options:
#   --sample_id        Output filename prefix (required)
#   --genome_version   Genome version used to exclude known problematic regions (default: GRCh38)
#
# Dependencies:
#   - bcftools must be available in the system path.
#   - A genome-specific exclusion list: resources/Genome_Regions_data.tsv
#     (should contain GRCh38 and/or GRCh37 region names to exclude from CNV inference)
#
# Example usage:
#   perl extract_baf_lrr.pl sample.gvcf.gz --sample_id sample_output --genome_version GRCh38
# ------------------------------------------


use Cwd 'abs_path';
use File::Basename;

my $script_dir  = dirname(abs_path($0)); # directory containing the script

# Command-line variable: --sample_id specifies the output file prefix
my $sample_id;
my $genome_version = "GRCh38";  # default value
my $GQ = "20";  # default value
my $DP = "10";  # default value
my $output_dir = ".";  # default: current directory

# Parse options: --sample_id and --genome_version
GetOptions(
    'sample_id=s'        => \$sample_id,
    'genome_version=s' => \$genome_version,
    'output_dir=s' => \$output_dir,
    'GQ=s' => \$GQ,
    'DP=s' => \$DP,
) or die "Usage: perl script.pl input.gvcf.gz --sample_id sample_id [--output_dir output_dir] [--genome_version GRCh37|GRCh38] [--GQ INT] [--DP INT]\n";

# Ensure one input file and --sample_id provided
@ARGV == 1 && $sample_id
  or die "Usage: perl script.pl input.gvcf.gz --sample_id sample_id [--output_dir output_dir] [--genome_version GRCh37|GRCh38]\n";

# Input gVCF file
my $input_gvcf = $ARGV[0];

# Global variable to store the mean coverage across SNPs
my $meancov;

# Step 1: Extract SNPs and calculate BAF
readVariantInfo("$output_dir/$sample_id.snp_metrics.tsv");

# Step 2: Compute Log R Ratio and generate final report
addLRRBAF("$output_dir/$sample_id.snp_metrics.tsv", "$output_dir/$sample_id.baf_lrr.tsv", $meancov, $sample_id, $output_dir);

# ===================================================
# FUNCTION: Extract SNPs and compute BAF
# Output: <prefix>.snp_metrics.tsv
# ===================================================
sub readVariantInfo {
    my ($readoutfile) = @_;

    # Build bcftools pipeline:
    my $command = "bcftools view -m3 -M3 -V indels $input_gvcf | " .               # Keep only biallelic SNPs
                  "bcftools filter -e 'format/GQ<$GQ|format/DP<$DP' | " .            # Filter out low-quality genotypes
                  "bcftools view -T ^$script_dir/resources/Genome_Regions_data_${genome_version}.tsv  | " .        # Exclude known problematic regions
                  "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO[\t%GT\t%AD\t%DP]\n' |";

    # Print the command for debugging
    print STDERR "DEBUG: Running command:\n$command\n";

    # Initialize counters
    my ($countsite, $countcov, $countsnp) = (0, 0, 0);

    # Open the bcftools pipeline
    open (VAR, $command) or die "Error: failed to execute bcftools pipeline: $!\n";

    # Open output file
    open (OUT, ">$readoutfile") or die "Error: cannot write to $readoutfile: $!\n";

    # Write header
    print OUT "Name\tCoverage\tBAC\tBAF\tLength\n";

    # Process each variant line
    while (<VAR>) {
        chomp;

        # Split the tab-delimited fields
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $GT, $AD, $DP) = split /\t/;

        # Skip entries with missing depth
        next if $DP eq ".";

        my $cov = $DP;
        $countsite++;
        $countcov += $cov;

        # Skip entries without AD field
        next if $AD eq ".";

        # AD format: ref_count, alt_count, nonref_count
        my ($AREF, $AALT, $NONREF) = split /,/, $AD;

        # Use alt_count as BAC (B Allele Count)
        my $bac = ($AALT ne ".") ? $AALT : next;

        # Compute B Allele Frequency
        my $baf = $bac / $cov;

        # Format position name: chr:start-end
        my $region = "$chr:$pos-$pos";

        # Write data line
        printf OUT "%s\t%d\t%d\t%.2f\t%d\n", $region, $cov, $bac, $baf, 1;

        $countsnp++;
    }

    # Compute mean coverage
    $meancov = $countcov / $countsite;

    print STDERR "NOTICE: Processed $countsite sites ($countsnp valid SNPs), mean coverage = $meancov\n";
}

# ===================================================
# FUNCTION: Compute Log R Ratio and write final output
# Output: <prefix>.baf_lrr.tsv (PennCNV and QuantiSNP compatible)
# ===================================================
sub addLRRBAF {
    my ($readinfile, $readoutfile, $meancov, $sample_id, $output_dir) = @_;

    # Fallback if mean coverage was not calculated
    $meancov ||= 30;

    print STDERR "NOTICE: Adding LRR/BAF to final report: $readoutfile\n";

    open (IN,  $readinfile)     or die "Error: cannot read $readinfile: $!\n";
    open (OUT, ">$readoutfile") or die "Error: cannot write $readoutfile: $!\n";

    # Read and validate header
    $_ = <IN>; chomp;
    /^Name\tCoverage/ or die "Error: invalid header in $readinfile: <$_>\n";

    # Write final header in PennCNV-compatible format
    print OUT "Name\tChr\tPosition\t$sample_id.Log R Ratio\t$sample_id.B Allele Freq\n";

    # Read each line of SNP data
    while (<IN>) {
        chomp;

        # Split input fields
        my ($name, $coverage, $bac, $baf, $length, $TYPE) = split /\t/;

        # Prevent division by zero
        $coverage = $coverage > 0 ? $coverage : 1;

        # Parse chromosome and position from name (e.g., chr1:12345-12345)
        my ($chr, $start, $end);
        if ($name =~ /^(?:chr)?(\w+):(\d+)-(\d+)$/) {
            ($chr, $start, $end) = ($1, $2, $3);
        } else {
            die "Error: invalid Name format <$name>\n";
        }

        # Calculate Log R Ratio = log(observed_coverage / mean_coverage)
        my $lrr = log($coverage / $meancov);

        # Write output line
        print OUT join("\t", $name, $chr, $start, $lrr, $baf), "\n";
    }
}


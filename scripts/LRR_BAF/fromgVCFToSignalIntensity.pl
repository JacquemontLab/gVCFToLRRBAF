#!/usr/bin/env perl
use warnings;     # Display helpful warnings for potential issues
use strict;       # Enforce variable declarations to avoid bugs
use Getopt::Long; # Module to parse command-line options (like --outfile)

# ------------------------------------------
# Script purpose:
# Given a compressed gVCF file (.g.vcf.gz), this script extracts biallelic SNPs,
# calculates B Allele Frequency (BAF) and Log R Ratio (LRR), and outputs two files:
# - <prefix>.read1       → Contains coverage, BAF, and B Allele Count (BAC)
# - <prefix>.finalReport → Format compatible with PennCNV and QuantiSNP
# ------------------------------------------

use Cwd 'abs_path';
use File::Basename;

my $script_dir  = dirname(abs_path($0)); # directory containing the script

# Command-line variable: --outfile specifies the output file prefix
my $outfile;
my $genome_version = "GRCh38";  # default value

# Parse options: --outfile and --genome_version
GetOptions(
    'outfile=s'        => \$outfile,
    'genome_version=s' => \$genome_version,
) or die "Usage: perl script.pl input.gvcf.gz --outfile outputprefix [--genome_version GRCh37|GRCh38]\n";

# Ensure one input file and --outfile provided
@ARGV == 1 && $outfile
  or die "Usage: perl script.pl input.gvcf.gz --outfile outputprefix [--genome_version GRCh37|GRCh38]\n";

# Input gVCF file
my $input_gvcf = $ARGV[0];

# Global variable to store the mean coverage across SNPs
my $meancov;

# Step 1: Extract SNPs and calculate BAF
readVariantInfo("$outfile.read1");

# Step 2: Compute Log R Ratio and generate final report
addLRRBAF("$outfile.read1", "$outfile.finalReport", $meancov);

# ===================================================
# FUNCTION: Extract SNPs and compute BAF
# Output: <prefix>.read1
# ===================================================
sub readVariantInfo {
    my ($readoutfile) = @_;

    # Build bcftools pipeline:
    my $command = "bcftools view -m3 -M3 -V indels $input_gvcf | " .               # Keep only biallelic SNPs
                  "bcftools filter -e 'format/GQ<20|format/DP<10' | " .            # Filter out low-quality genotypes
                  "bcftools view -T ^<(grep \"$genome_version\" $script_dir/resources/Genome_Regions_data.tsv)  | " .        # Exclude known problematic regions
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
# Output: <prefix>.finalReport (PennCNV-compatible)
# ===================================================
sub addLRRBAF {
    my ($readinfile, $readoutfile, $meancov) = @_;

    # Fallback if mean coverage was not calculated
    $meancov ||= 30;

    # Extract prefix from filename (e.g., "sample.read1" → "sample")
    my ($prefix) = $readinfile =~ m{([^/]+)\.read1$}
        or die "Error: cannot extract prefix from <$readinfile>\n";

    print STDERR "NOTICE: Adding LRR/BAF to final report: $readoutfile\n";

    open (IN,  $readinfile)     or die "Error: cannot read $readinfile: $!\n";
    open (OUT, ">$readoutfile") or die "Error: cannot write $readoutfile: $!\n";

    # Read and validate header
    $_ = <IN>; chomp;
    /^Name\tCoverage/ or die "Error: invalid header in $readinfile: <$_>\n";

    # Write final header in PennCNV-compatible format
    print OUT "Name\tChr\tPosition\t$prefix.Log R Ratio\t$prefix.B Allele Freq\n";

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


#!/usr/bin/env perl
use warnings;     # Display helpful warnings for potential issues
use strict;       # Enforce variable declarations to avoid bugs
use Getopt::Long; # Module to parse command-line options (like --sample_id)


# ------------------------------------------
# Script: fromgVCFToSignalIntensity.pl
#
# Purpose:
#   Processes a compressed gVCF file (.gvcf.gz) to:
#     1. Extract biallelic SNPs with adequate coverage and quality
#     2. Calculate per-SNP metrics: 
#        - Coverage (DP)
#        - B Allele Count (BAC)
#        - B Allele Frequency (BAF = BAC / DP)
#     3. Generate input for CNV detection tools (PennCNV, QuantiSNP):
#        - Log R Ratio (LRR) = log(observed coverage / mean coverage)
#        - BAF
#
# Output:
#   - <prefix>.snp_metrics.tsv : Coverage-based metrics per SNP
#   - <prefix>.baf_lrr.tsv     : Final table for CNV callers
#
# Required arguments:
#   <input.gvcf.gz>         Compressed gVCF input file (must be bgzipped and indexed)
#   --sample_id STRING      Output file prefix (e.g., "sample123")
#
# Optional arguments:
#   --genome_version STRING  Genome build: GRCh38 (default) or GRCh37
#   --output_dir STRING      Directory to save output files (default: current directory)
#   --GQ INT                 Minimum genotype quality (default: 20)
#   --DP INT                 Minimum depth of coverage (default: 10)
#
# Dependencies:
#   - bcftools must be installed and in the system PATH
#   - Exclusion list: resources/Genome_Regions_data_<genome_version>.tsv
#
# Example:
#   perl fromgVCFToSignalIntensity.pl sample.gvcf.gz \
#     --sample_id sample123 \
#     --genome_version GRCh38 \
#     --output_dir results/ \
#     --GQ 20 \
#     --DP 10 \
#     --bim .bimfile
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
my $bim_file;  #optional

# Parse options: --sample_id and --genome_version
GetOptions(
    'sample_id=s'        => \$sample_id,
    'genome_version=s' => \$genome_version,
    'output_dir=s' => \$output_dir,
    'GQ=s' => \$GQ,
    'DP=s' => \$DP,
    'bim=s' => \$bim_file,  # <-- new optional
) or die "Usage: perl script.pl input.gvcf.gz --sample_id sample_id [--output_dir output_dir] [--genome_version GRCh37|GRCh38] [--GQ INT] [--DP INT]\n";

# Ensure one input file and --sample_id provided
@ARGV == 1 && $sample_id
  or die "Usage: perl script.pl input.gvcf.gz --sample_id sample_id [--output_dir output_dir] [--genome_version GRCh37|GRCh38]\n";

# Input gVCF file
my $input_gvcf = $ARGV[0];

# ----- Build a set of (chr,pos) to keep if --bim is provided -----
my %keep;   # keys: "chrN\tPOS", e.g. "chr1\t10473"
if (defined $bim_file && length $bim_file) {
    open my $BIM, "<", $bim_file or die "Error: cannot read --bim $bim_file: $!\n";
    while (<$BIM>) {
        chomp;
        next if $_ eq "";
        my @f = split /\t/;          # BIM columns: CHR  ID  CM  POS ...
        next unless @f >= 4;
        my ($bchr, $bpos) = ($f[0], $f[3]);

        # normalize to WITH 'chr' (avoid double 'chr')
        $bchr =~ s/^chr//i;          # strip any existing 'chr' first
        $bchr = "chr$bchr";

        $keep{"$bchr\t$bpos"} = 1;
    }
    close $BIM;

}
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
                  "bcftools view -T ^$script_dir/resources/Genome_Regions_data_${genome_version}.bed  | " .        # Exclude known problematic regions
                  "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO[\t%GT\t%AD\t%DP]\n' |";

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
        
        # BIM whitelist: keep only positions present in the BIM set (WITH 'chr')
        if (%keep) {
            my $chr_norm = $chr;
            $chr_norm =~ s/^chr//i;      # enlève un éventuel 'chr' en tête
            $chr_norm = "chr$chr_norm";  # remet exactement un seul 'chr'
            next unless exists $keep{"$chr_norm\t$pos"};  # si absent du BIM -> on saute la ligne
        }

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
    $meancov = $countsite > 0 ? $countcov / $countsite : 30;

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


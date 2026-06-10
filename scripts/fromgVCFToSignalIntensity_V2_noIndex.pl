#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

# ------------------------------------------
# Script: fromgVCFToSignalIntensity.pl
#
# Purpose:
#   Processes a compressed gVCF file (.gvcf.gz) to extract BAF/LRR signals
#   for CNV detection with PennCNV or QuantiSNP.
#
#   Supports all gVCF formats:
#     - DeepVariant : all positions explicit (SNV + ref/ref)
#     - SNV-only    : only variant positions (e.g. All of Us)
#     - GATK gVCF   : reference blocks with END tag (e.g. HaplotypeCaller)
#
# Output:
#   - <prefix>.snp_metrics.tsv : coverage-based metrics per position
#   - <prefix>.baf_lrr.tsv     : BAF/LRR table (PennCNV/QuantiSNP compatible)
#
# Required arguments:
#   <input.gvcf.gz>        Compressed gVCF (bgzipped + indexed)
#   --sample_id STRING     Output file prefix
#
# Optional arguments:
#   --bim FILE             PLINK BIM file — positions to extract (recommended)
#   --genome_version STR   GRCh38 (default) or GRCh37
#   --output_dir STR       Output directory (default: current)
#   --GQ INT               Min genotype quality (default: 20)
#   --DP INT               Min depth of coverage (default: 10)
#
# Example:
#   perl fromgVCFToSignalIntensity.pl sample.gvcf.gz \
#     --sample_id SP0001234 --bim merged_dataset.bim \
#     --output_dir LRR_BAF/ --GQ 20 --DP 10
# ------------------------------------------

my $script_dir = dirname(abs_path($0));

my ($sample_id, $bim_file);
my $genome_version = "GRCh38";
my $min_GQ         = 20;
my $min_DP         = 10;
my $output_dir     = ".";

GetOptions(
    'sample_id=s'      => \$sample_id,
    'genome_version=s' => \$genome_version,
    'output_dir=s'     => \$output_dir,
    'GQ=i'             => \$min_GQ,
    'DP=i'             => \$min_DP,
    'bim=s'            => \$bim_file,
) or die "Usage: perl fromgVCFToSignalIntensity.pl input.gvcf.gz --sample_id ID [options]\n";

@ARGV == 1 && defined $sample_id
    or die "Usage: perl fromgVCFToSignalIntensity.pl input.gvcf.gz --sample_id ID [options]\n";

my $input_gvcf = $ARGV[0];

# ─── BIM data structures ──────────────────────────────────────────────────────
# %keep       : exact lookup   "chrN\tPOS" -> 1
# %chr_sorted : range queries  chrN -> [sorted positions]
my (%keep, %chr_sorted);

if (defined $bim_file && length $bim_file) {
    open my $BIM, "<", $bim_file or die "Cannot read BIM $bim_file: $!\n";
    while (<$BIM>) {
        chomp; next if /^\s*$/;
        my @f = split /\t/;
        next unless @f >= 4;
        my ($bchr, $bpos) = ($f[0], $f[3]);
        $bchr =~ s/^chr//i;
        $bchr = "chr$bchr";
        $keep{"$bchr\t$bpos"} = 1;
        push @{$chr_sorted{$bchr}}, $bpos;
    }
    close $BIM;
    $chr_sorted{$_} = [ sort { $a <=> $b } @{$chr_sorted{$_}} ]
        for keys %chr_sorted;
    my $n = scalar keys %keep;
    print STDERR "NOTICE: Loaded $n BIM positions\n";
}

my $meancov;
readVariantInfo("$output_dir/$sample_id.snp_metrics.tsv");
addLRRBAF(
    "$output_dir/$sample_id.snp_metrics.tsv",
    "$output_dir/$sample_id.baf_lrr.tsv",
    $meancov, $sample_id
);

# ─── Binary search: BIM positions in [start, end] ────────────────────────────
sub bim_in_range {
    my ($chr, $start, $end) = @_;
    my $arr = $chr_sorted{$chr} // [];
    return () unless @$arr;
    my ($lo, $hi) = (0, $#$arr);
    while ($lo < $hi) {
        my $mid = int(($lo + $hi) / 2);
        $arr->[$mid] < $start ? ($lo = $mid + 1) : ($hi = $mid);
    }
    my @result;
    while ($lo <= $#$arr && $arr->[$lo] <= $end) {
        push @result, $arr->[$lo++];
    }
    return @result;
}

# ─── Step 1: Extract positions and compute BAF ────────────────────────────────
sub readVariantInfo {
    my ($outfile) = @_;

    my $regions_bed = "$script_dir/resources/Genome_Regions_data_${genome_version}.bed";

    my $command =
        "bcftools view $input_gvcf | " .
        "bcftools filter -e 'format/GQ<$min_GQ|format/DP<$min_DP' | " .
        "bcftools view -T ^$regions_bed | " .
        "bcftools query -f " .
        "'%CHROM\t%POS\t%REF\t%ALT\t%INFO[\t%GT\t%AD\t%DP]\n' |";

    my ($countsite, $countcov, $countsnp) = (0, 0, 0);
    my $gatk_blocks = 0;

    open(VAR, $command) or die "bcftools pipeline failed: $!\n";
    open(OUT, ">$outfile") or die "Cannot write $outfile: $!\n";
    print OUT "Name\tCoverage\tBAC\tBAF\tLength\n";

    while (<VAR>) {
        chomp;
        my ($chr, $pos, $ref, $alt, $info, $GT, $AD, $DP) = split /\t/;

        # ── GATK reference block ─────────────────────────────────────────────
        if (($alt =~ /<NON_REF>|<\*>/) && ($info =~ /END=(\d+)/)) {
            my $block_end = $1;
            next if !defined $DP || $DP eq "." || $DP < $min_DP;

            if (!%keep) {
                warn "WARNING: GATK reference block at $chr:$pos-$block_end " .
                     "skipped (provide --bim to extract positions from blocks)\n"
                    unless $gatk_blocks++;
                next;
            }

            (my $chr_norm = $chr) =~ s/^chr//i;
            $chr_norm = "chr$chr_norm";

            for my $bpos (bim_in_range($chr_norm, $pos, $block_end)) {
                printf OUT "%s\t%d\t%d\t%.2f\t%d\n",
                    "$chr:$bpos-$bpos", $DP, 0, 0.0, 1;
                $countsite++;
                $countcov += $DP;
                $countsnp++;
            }
            next;
        }

        # ── Regular position: SNV or explicit ref/ref (DeepVariant / SNV-only) ─
        if (%keep) {
            (my $chr_norm = $chr) =~ s/^chr//i;
            $chr_norm = "chr$chr_norm";
            next unless exists $keep{"$chr_norm\t$pos"};
        }

        next if !defined $DP || $DP eq "." || $DP < $min_DP;

        my $cov = $DP;
        $countsite++;
        $countcov += $cov;

        my ($bac, $baf) = (0, 0);
        if (defined $AD && $AD ne ".") {
            my ($AREF, $AALT) = split /,/, $AD;
            $bac = (defined $AALT && $AALT ne ".") ? $AALT : 0;
            $baf = $cov > 0 ? $bac / $cov : 0;
        }

        printf OUT "%s\t%d\t%d\t%.2f\t%d\n",
            "$chr:$pos-$pos", $cov, $bac, $baf, 1;
        $countsnp++;
    }
    close VAR;
    close OUT;

    $meancov = $countsite > 0 ? $countcov / $countsite : 30;
    print STDERR "NOTICE: Processed $countsite sites ($countsnp valid SNPs), " .
                 "mean coverage = $meancov\n";
}

# ─── Step 2: Compute LRR and write final output ───────────────────────────────
sub addLRRBAF {
    my ($infile, $outfile, $meancov, $sample_id) = @_;
    $meancov ||= 30;

    print STDERR "NOTICE: Adding LRR/BAF to final report: $outfile\n";

    open(IN,  $infile)     or die "Cannot read $infile: $!\n";
    open(OUT, ">$outfile") or die "Cannot write $outfile: $!\n";

    $_ = <IN>; chomp;
    /^Name\tCoverage/ or die "Invalid header in $infile: <$_>\n";

    print OUT "Name\tChr\tPosition\t$sample_id.Log R Ratio\t$sample_id.B Allele Freq\n";

    while (<IN>) {
        chomp;
        my ($name, $coverage, $bac, $baf) = split /\t/;
        $coverage = $coverage > 0 ? $coverage : 1;

        $name =~ /^(?:chr)?(\w+):(\d+)-(\d+)$/
            or die "Invalid Name format: <$name>\n";
        my ($chr, $start) = ($1, $2);

        my $lrr = log($coverage / $meancov);
        print OUT join("\t", $name, $chr, $start, $lrr, $baf), "\n";
    }
    close IN;
    close OUT;
}

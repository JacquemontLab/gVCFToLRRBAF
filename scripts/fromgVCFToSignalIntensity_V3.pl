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
#     - GATK gVCF   : reference blocks with END tag (HaplotypeCaller)
#     - Any caller  : large reference blocks handled via binary search O(log n)
#
# Output:
#   - <prefix>.snp_metrics.tsv : coverage-based metrics per position
#   - <prefix>.baf_lrr.tsv     : BAF/LRR table (PennCNV/QuantiSNP compatible)
#
# Required arguments:
#   <input.gvcf.gz>        Compressed gVCF (bgzipped + tabix-indexed)
#   --sample_id STRING     Output file prefix
#
# Optional arguments:
#   --bim FILE             PLINK BIM file (recommended — positions to extract)
#   --genome_version STR   GRCh38 (default) or GRCh37
#   --output_dir STR       Output directory (default: current)
#   --GQ INT               Min genotype quality (default: 20)
#   --DP INT               Min depth of coverage (default: 10)
#
# Dependencies:
#   - bcftools >= 1.14 (for --regions-overlap)
#   - tabix-indexed gVCF (.tbi or .csi) when using --bim
#
# Example:
#   perl fromgVCFToSignalIntensity.pl sample.gvcf.gz \
#     --sample_id SP0001234 --bim merged_dataset.bim \
#     --output_dir LRR_BAF/ --GQ 20 --DP 10
# ------------------------------------------

my $script_dir = dirname(abs_path($0));

my ($sample_id, $bim_file);
my $genome_version = "GRCh38";
my $GQ_threshold   = 20;
my $DP_threshold   = 10;
my $output_dir     = ".";

GetOptions(
    'sample_id=s'      => \$sample_id,
    'genome_version=s' => \$genome_version,
    'output_dir=s'     => \$output_dir,
    'GQ=i'             => \$GQ_threshold,
    'DP=i'             => \$DP_threshold,
    'bim=s'            => \$bim_file,
) or die "Usage: perl fromgVCFToSignalIntensity.pl input.gvcf.gz --sample_id ID [options]\n";

@ARGV == 1 && defined $sample_id
    or die "Usage: perl fromgVCFToSignalIntensity.pl input.gvcf.gz --sample_id ID [options]\n";

my $input_gvcf = $ARGV[0];

# ─── BIM data structures ──────────────────────────────────────────────────────
# %keep       : exact lookup O(1)  →  "chrN\tPOS" -> 1
# %chr_sorted : range queries      →  chrN -> [sorted positions]
my (%keep, %chr_sorted);

my $targets_file = "$output_dir/$sample_id.bim_targets.tsv";
my $use_targets  = 0;

if (defined $bim_file && -s $bim_file) {
    print STDERR "NOTICE: Loading BIM file...\n";
    open my $BIM, "<", $bim_file or die "Cannot read BIM $bim_file: $!\n";
    open my $TGT, ">", $targets_file or die "Cannot write $targets_file: $!\n";

    while (<$BIM>) {
        chomp; next if /^\s*$/;
        my @f = split /\t/;
        next unless @f >= 4;
        my ($bchr, $bpos) = ($f[0], $f[3]);
        $bchr =~ s/^chr//i;
        $bchr = "chr$bchr";
        $keep{"$bchr\t$bpos"} = 1;
        push @{$chr_sorted{$bchr}}, $bpos;
        # 2 columns, 1-based — unambiguous for bcftools -R
        print $TGT "$bchr\t$bpos\n";
    }
    close $BIM;
    close $TGT;

    # Sort targets file (required by bcftools -R)
    system("sort -k1,1 -k2,2n -o $targets_file $targets_file") == 0
        or die "Failed to sort $targets_file\n";

    # Sort each chromosome's position array for binary search
    $chr_sorted{$_} = [ sort { $a <=> $b } @{$chr_sorted{$_}} ]
        for keys %chr_sorted;

    my $n = scalar keys %keep;
    print STDERR "NOTICE: Loaded $n BIM positions\n";
    $use_targets = 1;
}

my $meancov = 0;
readVariantInfo("$output_dir/$sample_id.snp_metrics.tsv");
addLRRBAF(
    "$output_dir/$sample_id.snp_metrics.tsv",
    "$output_dir/$sample_id.baf_lrr.tsv",
    $meancov, $sample_id
);

# ─── Binary search: BIM positions within [start, end] ────────────────────────
# O(log n + k) — efficient regardless of block size
sub bim_in_range {
    my ($chr, $start, $end) = @_;
    my $arr = $chr_sorted{$chr} // [];
    return () unless @$arr;

    # Find first index >= $start
    my ($lo, $hi) = (0, $#$arr);
    while ($lo < $hi) {
        my $mid = int(($lo + $hi) / 2);
        $arr->[$mid] < $start ? ($lo = $mid + 1) : ($hi = $mid);
    }
    return () if $arr->[$lo] < $start;
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

    # Separate pipes: -R and -T ^ must not be mixed in the same bcftools call
    # --regions-overlap record: returns gVCF blocks even when POS < query position
    # Requires bcftools >= 1.14 and tabix-indexed gVCF
    my $view_cmd = $use_targets
        ? "bcftools view -R $targets_file --regions-overlap record $input_gvcf"
        : "bcftools view $input_gvcf";

    my $command =
        "$view_cmd | " .
        "bcftools view -T ^$regions_bed | " .
        "bcftools query -f '%CHROM\t%POS\t%INFO/END\t%REF\t%ALT\t[%GQ\t%DP\t%MIN_DP\t%AD]\n' |";

    print STDERR "NOTICE: Running bcftools extraction...\n";
    open(VAR, $command) or die "bcftools pipeline failed: $!\n";
    open(OUT, ">$outfile") or die "Cannot write $outfile: $!\n";
    print OUT "Name\tCoverage\tBAC\tBAF\tLength\n";

    my ($countsnp, $sumcov) = (0, 0);

    while (<VAR>) {
        chomp;
        my ($chr, $pos, $end, $ref, $alt, $gq, $dp, $min_dp, $ad) = split /\t/;

        # Normalize chr prefix
        (my $chr_norm = $chr) =~ s/^chr//i;
        $chr_norm = "chr$chr_norm";

        # Resolve block end (or single position if no END tag)
        $end = ($end ne "." && $end =~ /^\d+$/) ? $end : $pos;

        # ── Determine positions to process ───────────────────────────────────
        my @positions_to_process;

        if ($pos == $end) {
            # Single position (SNV or single-base ref/ref): O(1) lookup
            if (!$use_targets || exists $keep{"$chr_norm\t$pos"}) {
                @positions_to_process = ($pos);
            }
        } else {
            # Reference block: binary search — O(log n + k) for any block size
            if ($use_targets) {
                @positions_to_process = bim_in_range($chr_norm, $pos, $end);
            } else {
                # No BIM: unavoidable linear expansion
                warn "WARNING: large reference block $chr_norm:$pos-$end without --bim\n"
                    if ($end - $pos) > 10000;
                @positions_to_process = ($pos .. $end);
            }
        }

        next unless @positions_to_process;

        # ── Resolve coverage: GATK (DP) vs DeepVariant/other (MIN_DP) ────────
        my $cov = 0;
        if (defined $dp && $dp ne "." && $dp > 0) {
            $cov = $dp;
        } elsif (defined $min_dp && $min_dp ne "." && $min_dp > 0) {
            $cov = $min_dp;
        } else {
            next;
        }

        # ── Quality filters in Perl (after DP/MIN_DP resolution) ─────────────
        my $qual = (defined $gq && $gq ne ".") ? $gq : 0;
        next if $qual < $GQ_threshold;
        next if $cov  < $DP_threshold;

        # ── BAC: sum ALL alt alleles (bi-allelic model for PennCNV) ──────────
        my ($bac, $baf) = (0, 0);
        if (defined $ad && $ad ne ".") {
            my @alleles = split /,/, $ad;
            for my $i (1 .. $#alleles) {
                $bac += $alleles[$i]
                    if defined $alleles[$i] && $alleles[$i] ne ".";
            }
            $baf = $cov > 0 ? $bac / $cov : 0;
            $baf = 1 if $baf > 1;
        }

        # ── Write one line per BIM position ──────────────────────────────────
        for my $current_pos (@positions_to_process) {
            printf OUT "%s\t%d\t%d\t%.4f\t%d\n",
                "$chr_norm:$current_pos-$current_pos", $cov, $bac, $baf, 1;
            $countsnp++;
            $sumcov += $cov;
        }
    }
    close VAR;
    close OUT;

    # Mean coverage computed on exported SNPs only (not on raw gVCF lines)
    $meancov = $countsnp > 0 ? $sumcov / $countsnp : 30;
    print STDERR "NOTICE: Processed $countsnp valid SNPs, mean coverage = $meancov\n";

    unlink $targets_file if -e $targets_file;
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

        print OUT join("\t", $name, $chr, $start, log($coverage / $meancov), $baf), "\n";
    }
    close IN;
    close OUT;
}

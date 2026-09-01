# gVCFToLRRBAF

A Nextflow DSL2 pipeline for generating Log R Ratio (LRR) and B Allele Frequency (BAF) signals from whole-genome sequencing (WGS) gVCF files for downstream copy number variant (CNV) detection with established LRR/BAF-based callers.

> **A preprint describing this pipeline has been submitted to bioRxiv.**

---

## Overview

gVCFToLRRBAF takes population-level VCF (pVCF) files and individual-level gVCF files as input and generates per-sample LRR and BAF signals at selected SNV positions.

Population-level VCF files are first used to construct a cohort-level SNV reference panel after genotype quality-control filtering. The resulting PLINK BIM file defines the SNV positions retained for downstream signal extraction.

For each individual, LRR and BAF signals are then derived from the gVCF at the retained SNV positions. A subset of unrelated individuals is selected, and their generated BAF values are used to estimate a cohort-specific Population Frequency of B Allele (PFB) at each retained SNV position.

The generated LRR and BAF signals can subsequently be used for CNV calling with established LRR/BAF-based callers. QuantiSNP uses the generated LRR/BAF signals, whereas PennCNV uses both the LRR/BAF signals and the cohort-specific PFB.

---

## Workflow overview

```text
Population VCFs
      │
      ▼
PLINK conversion and QC
      │
      ▼
Cohort-level SNV reference panel
      │
      ▼
Final BIM file
      │
      ▼
Individual gVCFs
      │
      ▼
LRR / BAF signal extraction
      │
      ├──────────────────────────────► QuantiSNP
      │
      ├──────────────────────────────► LRR / BAF signals ─────────────┐
      │                                                               │
      ▼                                                               │
Selected unrelated individuals                                       │
      │                                                               │
      ▼                                                               │
BAF aggregation across samples                                       │
      │                                                               │
      ▼                                                               │
Cohort-specific PFB ──────────────────────────────────────────────────┤
                                                                      ▼
                                                                   PennCNV
```

For QuantiSNP, CNV calling is performed using the generated LRR and BAF signals. For PennCNV, CNV calling uses the generated LRR and BAF signals together with the cohort-specific PFB.

---

## Key features

- Generates LRR and BAF signals directly from WGS gVCF files
- Uses a cohort-level SNV reference panel derived from population-level VCF data
- Restricts signal extraction to retained biallelic SNV positions
- Estimates a cohort-specific PFB from generated BAF signals in selected unrelated individuals
- Produces signals compatible with PennCNV and QuantiSNP
- Tested with DeepVariant and GATK HaplotypeCaller gVCFs, as well as SNV-only VCF representations
- Implemented in Nextflow DSL2
- Supports parallel execution on SLURM HPC clusters
- Processes large cohorts using batched execution
- Includes automatic retry with increased resource allocation for selected processes

---

## Requirements

### Software

| Tool | Version tested / required | Purpose |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.10 | Workflow management |
| [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/) | 2.0 | Genotype filtering and merging |
| [bcftools](https://samtools.github.io/bcftools/) | ≥ 1.14 | VCF/gVCF processing |
| [Perl](https://www.perl.org/) | ≥ 5.30 | LRR/BAF signal extraction |
| [GNU Parallel](https://www.gnu.org/software/parallel/) | ≥ 20230522 | Parallel sample processing |
| [Python](https://www.python.org/) | ≥ 3.10 | PFB computation |
| [Polars](https://pola.rs/) | ≥ 1.0 | Dataframe operations |

---

## Input files

### Population-level VCF files

One population-level VCF file per autosome is expected for construction of the cohort-level SNV reference panel.

Example:

```text
cohort.chr1.vcf.gz
cohort.chr2.vcf.gz
...
cohort.chr22.vcf.gz
```

The files should contain cohort genotype information and be compatible with PLINK 2.0.

### Individual gVCF files

One compressed gVCF file per individual is required for LRR/BAF signal extraction.

Example:

```text
SAMPLE001.g.vcf.gz
SAMPLE002.g.vcf.gz
```

The workflow has been tested with:

- DeepVariant gVCFs
- GATK HaplotypeCaller gVCFs
- SNV-only VCF representations

---

## Cohort-level SNV reference panel

Population-level VCF files from chromosomes 1–22 are converted to PLINK format and merged.

The genotype quality-control filters used in the manuscript are:

```text
--geno 0.05
--mind 0.05
--hwe 1e-6
--maf 0.001
```

These correspond to:

- variant missingness ≤ 5%
- individual missingness ≤ 5%
- Hardy-Weinberg equilibrium p ≥ 1 × 10⁻⁶
- minor allele frequency ≥ 0.001

Biallelic SNVs passing these filters are retained.

Variants overlapping problematic genomic regions are excluded before construction of the final cohort-level SNV reference panel.

The resulting PLINK BIM file contains the retained SNV positions used for downstream LRR/BAF signal extraction.

---

## Problematic genomic regions

Problematic genomic regions used for filtering include:

- segmental duplications
- centromeres
- telomeres
- the major histocompatibility complex (MHC)
- UCSC problematic-region tracks, including UCSC Unusual Regions, ENCODE Blacklist v2, and GRC exclusion regions

Coordinates were retrieved from the UCSC Genome Browser and Genome Reference Consortium resources in April 2025. Overlapping intervals were merged using BEDTools.

The annotation files, source URLs, and commands used to generate these resources are available in:

```text
scripts/resources/
```

Resources are provided for both GRCh37 and GRCh38.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/zedicush/gVCFToLRRBAF.git
cd gVCFToLRRBAF
```

Make sure Nextflow and the required dependencies are available in your environment.

---

## Configuration

Pipeline parameters can be defined in `nextflow.config` and overridden at runtime.

Example:

```groovy
params {

    inputDirpVCF = "/path/to/pvcf_directory"
    inputDirgVCF = "/path/to/gvcf_directory"
    outputDir    = "/path/to/output_directory"

    topN              = 1000
    gvcf_batch_size   = 125
    filter_batch_size = 250
    n_batches         = 20
    polars_threads    = 4

    PLINK_filter_options =
        "--geno 0.05 --mind 0.05 --hwe 1e-6 --maf 0.001"
}
```

### Main parameters

| Parameter | Description |
|---|---|
| `inputDirpVCF` | Directory containing population-level VCF files |
| `inputDirgVCF` | Directory containing individual gVCF files |
| `outputDir` | Output directory |
| `topN` | Number of unrelated individuals used for PFB estimation |
| `gvcf_batch_size` | Number of gVCFs processed per batch |
| `filter_batch_size` | Number of LRR/BAF files processed per filtering batch |
| `n_batches` | Number of batches used for PFB computation |
| `polars_threads` | Number of threads used by Polars |
| `PLINK_filter_options` | PLINK quality-control parameters |

---

## SLURM configuration

The workflow supports parallel execution on SLURM-based HPC systems.

Set the appropriate SLURM account in `nextflow.config`:

```groovy
process {
    clusterOptions = '--account=YOUR_ACCOUNT'
}
```

Resource requirements can be adapted to the target HPC environment.

---

## Usage

### Run the workflow

```bash
nextflow run main.nf -with-report
```

### Run through SLURM

Example:

```bash
sbatch \
    --account=YOUR_ACCOUNT \
    --time=48:00:00 \
    --mem=4G \
    --cpus-per-task=2 \
    --wrap="module load nextflow && nextflow run main.nf -with-report"
```

### Resume an interrupted run

Nextflow caching allows interrupted workflows to be resumed:

```bash
nextflow run main.nf -resume -with-report
```

### Override parameters at runtime

Example:

```bash
nextflow run main.nf \
    --inputDirgVCF /path/to/gvcf \
    --inputDirpVCF /path/to/pvcf \
    --outputDir results \
    --topN 500 \
    --gvcf_batch_size 100
```

---

## Pipeline steps

| Step | Process | Description |
|---|---|---|
| 1 | `pVCF_to_plink` | Convert population-level VCF files from chromosomes 1–22 to PLINK format |
| 2 | `plink_filter_merge` | Apply genotype QC filters and merge autosomes |
| 3 | `select_unrel_topN` | Select unrelated individuals for PFB estimation |
| 4 | `gvcf_to_signalintensity` | Extract LRR and BAF signals from individual gVCF files |
| 5 | `split_iids` | Divide selected individuals into batches |
| 6 | `pfb_partial` | Compute partial PFB statistics from generated BAF signals |
| 7 | `pfb_reduce` | Aggregate partial estimates into the final PFB file |
| 8 | `filter_baf_lrr_files` | Filter LRR/BAF files using the retained cohort-level reference positions |

---

## LRR and BAF signal extraction

For each individual, signal extraction is restricted to biallelic SNV positions included in the cohort-level SNV reference panel.

Information contained in the gVCF, including genotype, total depth, and allele depth, is used to derive LRR and BAF metrics.

The primary signal extraction script is:

```text
scripts/fromgVCFToSignalIntensity.pl
```

The workflow generates one LRR/BAF signal file per individual.

---

## Population Frequency of B Allele (PFB)

A cohort-specific Population Frequency of B Allele (PFB) file is estimated from the BAF signals generated by the workflow.

A subset of unrelated individuals is selected from the cohort, and their per-sample BAF values are aggregated at each SNV position retained in the cohort-level reference panel.

PFB estimation is performed in parallel using Python and Polars:

```text
scripts/pfb_partial.py
scripts/pfb_reduce.py
```

Partial estimates are calculated across batches of selected individuals and subsequently combined to generate:

```text
PFB.tsv
```

The resulting PFB represents the cohort-level frequency of the B allele at each retained SNV position.

The PFB file is used as an additional input for PennCNV and is not required for QuantiSNP.

---

## Output

Example output structure:

```text
results/
├── merged_dataset.bed
├── merged_dataset.bim
├── merged_dataset.fam
│
├── stats/
│   ├── merged_dataset.missing.smiss
│   └── merged_dataset.missing.vmiss
│
├── keep_unrel_topN.iid.txt
│
├── LRR_BAF/
│   └── <SAMPLE_ID>.baf_lrr.tsv
│
├── PFB.tsv
│
└── LRR_BAF_filtered/
    └── <SAMPLE_ID>_filtered.baf_lrr.tsv
```

### Main outputs

| Output | Description |
|---|---|
| `merged_dataset.bim` | Cohort-level SNV reference panel |
| `keep_unrel_topN.iid.txt` | Individuals selected for PFB estimation |
| `LRR_BAF/<SAMPLE_ID>.baf_lrr.tsv` | Per-sample LRR and BAF signals |
| `PFB.tsv` | Cohort-specific Population Frequency of B Allele file for PennCNV |
| `LRR_BAF_filtered/<SAMPLE_ID>_filtered.baf_lrr.tsv` | Filtered per-sample LRR/BAF signals |

---

## gVCF format compatibility

The signal extraction workflow has been tested with:

| Input representation | Description |
|---|---|
| DeepVariant gVCF | gVCF generated using DeepVariant |
| GATK HaplotypeCaller gVCF | gVCF containing reference blocks |
| SNV-only VCF | Variant-only representation without reference blocks |

Because gVCF structures and FORMAT fields can differ between variant-calling pipelines, users should verify compatibility when applying the workflow to other gVCF implementations.

---

## Downstream CNV calling

gVCFToLRRBAF generates LRR and BAF signals for downstream CNV calling with established LRR/BAF-based callers.

The manuscript evaluates:

- [PennCNV](https://penncnv.openbioinformatics.org/)
- QuantiSNP

For PennCNV, CNV calling uses the generated per-sample LRR/BAF signals together with the cohort-specific PFB.

QuantiSNP uses the generated LRR/BAF signals and does not require the PFB file.

Installation, configuration, and licensing of PennCNV and QuantiSNP should follow the documentation provided by their respective developers.

---

## Computational scalability

The workflow was evaluated on large WGS cohorts.

In the SPARK cohort, LRR/BAF signals were generated for 12,509 individuals at an average of approximately 2.7 million SNV positions per individual. Signal extraction was completed in approximately 4 hours using 192 CPUs.

The workflow was subsequently applied to 414,000 individuals from the All of Us cohort. The analysis was completed in approximately 96 hours of wall-clock time.

Resource requirements depend on cohort size, gVCF representation, batch size, storage performance, and HPC configuration.

---

## Citation

A preprint describing gVCFToLRRBAF has been submitted to bioRxiv.

The full citation and DOI will be added once the preprint is publicly available.

---

## License

A software license will be added before the manuscript release.

Please refer to the `LICENSE` file once the final license has been selected.

---

## Contact

**Mame Seynabou Diop**  
Université de Montréal / CHU Sainte-Justine Research Centre  
Email: `mame.seynabou.diop@umontreal.ca`

**Jacquemont Lab**  
CHU Sainte-Justine Research Centre  
Montréal, Québec, Canada

# Perl Script fromgVCFToSignalIntensity.pl: Extracting BAF and Log R Ratio from a gVCF File

## Description

This Perl script,fromgVCFToSignalIntensity.pl, extracts SNP-level information from a gVCF file, including coverage, B Allele Frequency (BAF), and Log R Ratio (LRR).

It relies on bcftools to filter variants, exclude problematic regions, and reformat the output into a structure compatible with CNV analysis tools.

## Output Files

- <prefix>.read1 : Contains SNPs with coverage (Coverage), B Allele Count (BAC), and BAF.
- <prefix>.finalReport : Contains the same SNPs with the following columns: **Name, Chr, Position, Log R Ratio** (calculated from global mean coverage), and **B Allele Frequency**. 

This format is **compatible with PennCNV and QuantiSNP**, allowing direct downstream integration.

## Requirements

- Perl
- bcftools installed and accessible in your system PATH
- A compressed and indexed gVCF file (.gvcf.gz)

## Usage
perl fromgVCFToSignalIntensity.pl sample_id.gvcf.gz --outfile sample_id

Example Output (sample_id.finalReport)

Name               Chr   Position   sample_id.Log R Ratio   sample_id.B Allele Freq
chr1:10000-10000    1     10000              -0.08                    0.50
chr1:10234-10234    1     10234               0.12                    0.33


Processing Steps

SNP Filtering and Extraction:
	•	Keep only biallelic variants
	•	Exclude INDELs
	•	Keep variants with DP >=10 and GQ >=20
	•	Exclude known problematic regions (e.g., centromeres, segmental duplications, HLA)
 
BAF Calculation:
	•	BAF = alt allele count / total depth
 
Log R Ratio Calculation:
	•	LRR = log(observed coverage / mean coverage )


## Integration into a Parallel Pipeline (LRR/BAF Extraction Pipeline (MPI-Parallelized))

This script is part of a larger pipeline designed to extract signal intensity (BAF and LRR) from multiple individuals **in parallel** using **MPI** (OpenMPI + mpi4py) and the SLURM job scheduler.

The pipeline consists of the following components:

1. **`launchAllFinalReport.sh`**  
   SLURM job script that:
   - Sets up the working environment
   - Requests resources
   - Launches the pipeline with `mpiexec` on `mpiFinalReport.py`

2. **`mpiFinalReport.py`**  
   Python script that:
   - Splits gVCF files across MPI ranks (one sample per process)
   - Launches `launchFinalReport.sh` for each sample
   - Ensures efficient parallelism across nodes

3. **`launchFinalReport.sh`**  
   Bash script that:
   - Takes a sample name and input gVCF path
   - Sets output file names and paths
   - Calls the Perl script `fromgVCFToSignalIntensity.pl` with proper arguments

4. **`fromgVCFToSignalIntensity.pl`**  
   Perl script that:
   - Filters **biallelic SNPs** (`DP ≥ 10`, `GQ ≥ 20`)
   - Excludes **problematic regions**:
     - Centromeres
     - Segmental duplications
     - HLA
   - Computes:
     - **BAF (B Allele Frequency)**
     - **LRR (Log R Ratio)**
   - Outputs:
     - `.read1` file (not used for CNV analysis) and  `.finalReport` file (LRR and BAF for each SNV)

##  How to Run 
1. Prepare a **sample list** (one sample ID per line).
2. Make sure your input gVCFs follow a standard naming format (e.g., `<sample_id>.g.vcf.gz`).
3. Submit the MPI pipeline via SLURM:

sbatch launchAllFinalReport.sh

## output
For each sample, you will obtain:
	•	<sample_id>.read1 
	•	<sample_id>.finalReport 

## Dependencies
Perl (v5 or higher)
Python 3 with mpi4py installed
OpenMPI
SLURM scheduler

 
Author
Mame Seynabou Diop PhD Candidate in Bioinformatics 2025

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
Author
Mame Seynabou Diop PhD Candidate in Bioinformatics 2025

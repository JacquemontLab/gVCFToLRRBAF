# PLINK Preparation Module – WGS CNV Pipeline

This module handles the conversion of pVCF (Population VCF) files to PLINK binary format and applies quality control filters.

## 1. Convert pVCF to PLINK (pVCF_to_plink.sh)

This script converts a pVCF file into PLINK binary files (.bed, .bim, .fam).  
Each pVCF file contains all variants for one chromosome (here across 12,509 individuals).

### File organization

- **Input directory (inputDir)**:  
  Contains 24 files, one per chromosome.  
  Example: wgs_12509_genome.deepvariant.chr1.vcf.gz

- **Output directory (outputDir)**:  
  $workDir/pvcfToPlink 
  Example outputs:  
  - wgs_12509_genome.deepvariant.chr1.bed 
  - wgs_12509_genome.deepvariant.chr1.bim  
  - wgs_12509_genome.deepvariant.chr1.fam

### Required modules

- StdEnv/2020 (Compute Canada environment)
- plink/2.00a3.6 – [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/)

### Run the script

To convert one chromosome:
sbatch pVCF_to_plink.sh <CHR>
Example: sbatch pVCF_to_plink.sh 1


## 2. Process all chromosomes (launch_allpVCF.sh) at once
If you want to process all chromosomes (1-22, X, Y) at once, you can use the launch_allpVCF.sh script:
./launch_allpVCF.sh

This will submit one job per chromosome using sbatch pVCF_to_plink.sh <CHR>

## 3. Apply PLINK filters (plink_filters.sh)
Once the PLINK files are generated, quality control filters are applied:
	•	--geno 0.05 : remove SNVs with >5% missingness
	•	--mind 0.05 : remove individuals with >5% missing genotypes
	•	--hwe 0.000001 : remove SNVs out of Hardy-Weinberg equilibrium
	•	--maf 0.0001 : remove SNVs with minor allele frequency < 1/1000
Each chromosome is filtered independently.

## 4. Merge and Missingness Analysis
After filtering:
	•	Merge all filtered chromosomes using plink2 --pmerge-list
	•	Compute missing data statistics

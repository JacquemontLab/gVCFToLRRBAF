pVCF_to_plink.sh : Convert a pVCF (Population VCF) to PLINK binary format (BED/BIM/FAM)

This script converts a pVCF file (Population VCF) into PLINK binary files (.bed, .bim, .fam).

The script is designed to be launched independently for each chromosome to allow better resource management on a computing cluster.

Indeed, pVCF files are very large because they contain all variants for an entire chromosome across a large number of individuals (in this case, 12,509). To efficiently manage memory and processing time, each chromosome is processed separately.

File organization

     Input directory (inputDir): /home/$USER/projects/rrg-jacquese/All_user_common_folder/RAW_DATA_old/SPARK_V2_V3-2021_10/iWGS_v1_gVCF/deepvariant/pvcf This directory contains one pVCF file per chromosome (24 files total). Example: wgs_12509_genome.deepvariant.chr1.vcf.gz
     Output directory (outputDir): $workDir/pvcfToPlink This directory will contain the resulting PLINK files. Example: wgs_12509_genome.deepvariant.chr1.bed wgs_12509_genome.deepvariant.chr1.bim wgs_12509_genome.deepvariant.chr1.fam

Required module
StdEnv/2020. # compute canada environmmenet 2020
plink/2.00a3.6 ### PLINK 2 (https://www.cog-genomics.org/plink/2.0/)


How to run the script
sbatch pVCF_to_plink.sh <CHR>

Example: to process chromosome 1
sbatch pVCF_to_plink.sh 1

Process all chromosomes
If you want to process all chromosomes (1-22, X, Y) at once, you can use the launch_allpVCF.sh script:
./launch_allpVCF.sh

This will submit one job per chromosome using sbatch pVCF_to_plink.sh <CHR>

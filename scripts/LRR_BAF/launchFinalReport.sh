#!/bin/bash

# ----------- INPUT -----------
inputSample=$1

# ----------- DIRECTORIES -----------
workDir=$PWD

outputDir=$workDir/finalReport
mkdir -p "$outputDir"


# ----------- SAMPLE ID EXTRACTION -----------
sample_base=$(basename "$inputSample")
sampleId=${sample_base%.gvcf*}

# ----------- RUN THE CONVERSION from gVCF to LRR and BAF-----------
perl $workDir/fromgVCFToSignalIntensity.pl "$inputSample" --outfile "$outputDir/$sampleId"

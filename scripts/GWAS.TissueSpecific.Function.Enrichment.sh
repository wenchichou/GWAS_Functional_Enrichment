#!/bin/bash
# Tissue-Specific Functional Enrichment Analysis for GWAS hits
# This pipeline use HaploReg's ENCODE (ChromHMM) annotation to performe hypergeometric test for each SNP using its SNPs in LC block

use Python-2.7
# required python packages: bs4 

# enviroment setup
#resultDIR=$1 # where you store the intermediate results 
resultDIR="/home/unix/wcchou/gsapWenChi/tmp"
#inputSNPList=$2 # the list of input SNPs
inputSNPList="/home/unix/wcchou/gsapWenChi/tmp/SNPs_grepStrength.txt"
#collectSNPinLD_Py=$3 # a web parser to get SNPs in LD block by a given SNP
collectSNPinLD_Py="/home/unix/wcchou/gsapWenChi/tmp/collectSNPinLD.py"
#collectChromHMM_Sh=$4 # a web parser to collect ChromHMM annotations
collectChromHMM_Sh="/home/unix/wcchou/gsapWenChi/tmp/collectChromHMM.sh"
#enrichmentTest_R=$5 # an R script to run enrichmet analysis
enrichmentTest_R="/home/unix/wcchou/gsapWenChi/tmp/enrichmentTest.TissueSpecific.GWAS.r"
#enrichmentTestFun_R=$6 # an R script contains the functions for enrichment analysis
enrichmentTestFun_R="/home/unix/wcchou/gsapWenChi/tmp/enrichmentTest.TissueSpecific.GWAS.function.r"
#cellTable=$7 # a table describing the categories of all cell types
cellTable="/home/unix/wcchou/gsapWenChi/tmp/allCell.group.ID.description.txt"

# change dir to the result directory 
cd ${resultDIR} 

# 1) create folders for each input SNP
# 2) collect SNPs in the LD block of the given SNP [a python parser]
# 3) collect ChromHMM annotations of all SNPs [a shell parser]
for SNP in `cat ${inputSNPList}`; do 
		echo "working on ${SNP}";
		mkdir ${SNP}; 
		cd ${SNP};
		python ${collectSNPinLD_Py} ${SNP} SNPsInLDBlock.txt 
		bash ${collectChromHMM_Sh} SNPsInLDBlock.txt
		cd ${resultDIR};
done

# 4) perform hypergeometric test to find enriched regulatory functions in the given cell types [a R script]
# Rscript --vanilla Rscript RscripFunction DIR_collectedChromHMM enrichmentTable inputSNP iterationForPermutation cellTable
Rscript --vanilla ${enrichmentTest_R} ${enrichmentTestFun_R} ${resultDIR} gwasEnrichmentResults.txt ${inputSNPList} 100000 ${cellTable} 


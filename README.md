# GWAS_Functional_Enrichment
This is a pipeline designed for Tissue-specific regulatory-element enrichment analyses of the GWAS loci.
It will generate hypergeometric test p-values for all tissue-specificity by a set of given SNPs

The pipeline was first used at test GWAS loci associated with lean body mass. The work was first published on Nature Communication 2017; 8: 80. (Large meta-analysis of genome-wide association studies identifies five loci for lean body mass)
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5517526/table/Tab2/

The pipeline was also applied on other projects checking 
1) GWAS loci associated with handgrip and lower body strength (Aging Cell. 2016 Oct; 15(5): 792–800.)
2) new GWAS loci associated with lean body mass (submitted to Nature Communication)

Usage:
1) Clone this repository
2) install bs4, a python package
3) edit the environment setting in the schell script named "GWAS.TissueSpecific.Function.Enrichment.sh" in the "scripts" directory

3.1) set filepath of the input SNP list

3.2) set filepath of the result folder

4) run

(run the default setting and find your results at ./results/enrichment_grepStrenth/gwasEnrichmentResults.txt)

bash ./scripts/GWAS.TissueSpecific.Function.Enrichment.sh


posted by Wen-Chi Chou
#1/31/2018

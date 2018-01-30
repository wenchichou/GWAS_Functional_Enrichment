#!/bin/bash
# Wen-Chi Chou Nov. 19, 2013
# This SHELL script uses a list of SNPs (rs number) in the same LD block to get ENCODE/ChromHMM annotations from haploreg

# Test SNPlistFilePath="/home/wenchichou/project/leanmass/docs/leanMass_commonVariants/rs2943656/snpList"
# pass the path of the SNP list file. Format: one rs number per line
SNPlistFilePath=$1

# web parser
# get ChromHMM annotations from HaploReg
for i in `cat $SNPlistFilePath`; do 
	url=`echo "http://archive.broadinstitute.org/mammals/haploreg/detail_v4.php?query=&id="$i`
	echo $i;
	# search Haploreg to get ChromHMM annotations
	curl $url| sed -e 's/<\/table><p>/\n/g' -e 's/Regulatory /\nRegulatory /g'| grep "Regulatory chromatin states"| sed -e 's/^/\t/g' -e 's/<tr><td>/\t/g' -e 's/<\/b><p><table span class="detailtable">\t<b>/\t/g' -e 's/<\/td><td><b>/\t/g' -e 's/<\/td>$//g' -e 's/<\/td><td>/\t/g' -e 's/<\/td><\/tr>/\n/g'| grep -v "Regulatory"| grep -v "^$"| sed -e 's/^\t//g'| grep -v "</table>"| sed -e 's/<\/td><td /\t/'| grep 'Enh\|_Pro\|DNase'| cut -f1-4 > $SNPlistFilePath.$i.haploreg4.ChromHMM 
	sleep 1 # DO NOT CRASH THE SERVER
done

mkdir -p ./ChromHMM_DIR/
# move all ChromHMM tables to the same directory
mv ./*ChromHMM ./ChromHMM_DIR/
# collect all cell types from all ChromHMM tables
# the cell list will be used in the hypergeometric test 
cat ./ChromHMM_DIR/*ChromHMM| cut -f4 > ./ChromHMM.Cell


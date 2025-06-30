#!/bin/bash

# Names can't be more than 15 characters long
# $1 = tab-delimited list. Name first column, aa sequence second column
# $2 = path to the output directory

cat $1 | \
while read line
do
  output=`echo $line | awk '{print $1}'`
  seq=`echo $line | awk '{print $2}'`
	path=$2

	echo "solutions 1" > ${path}/${output}.dnaworks.job
	echo "repeat 8" >> ${path}/${output}.dnaworks.job
	echo "LOGFILE \"${path}/${output}.dnaworks.out\"" >> ${path}/${output}.dnaworks.job
	echo "pattern" >> ${path}/${output}.dnaworks.job
	echo "  BamHI GGATCC" >> ${path}/${output}.dnaworks.job
	echo "  NdeI CATATG" >> ${path}/${output}.dnaworks.job # NdeI
	echo "  XhoI CTCGAG" >> ${path}/${output}.dnaworks.job # XhoI
	echo "  NheI GCTAGC" >> ${path}/${output}.dnaworks.job # NheI
	echo "  BsaI GGTCTC" >> ${path}/${output}.dnaworks.job # BsaI
	echo "  BsaI GAGACC" >> ${path}/${output}.dnaworks.job # BsaI
	echo "  PolyA AAAAAAAA" >> ${path}/${output}.dnaworks.job
	echo "  PolyG GGGGG" >> ${path}/${output}.dnaworks.job
	echo "  PolyT TTTTTTTT" >> ${path}/${output}.dnaworks.job
	echo "  PolyC CCCCCCCC" >> ${path}/${output}.dnaworks.job
	echo "  Aarl  CACCTGC" >> ${path}/${output}.dnaworks.job
	echo "//" >> ${path}/${output}.dnaworks.job
	echo >> ${path}/${output}.dnaworks.job
  echo "codon codon S. cerevisiae" >> ${path}/${output}.dnaworks.job
	echo "protein" >> ${path}/${output}.dnaworks.job
	echo $seq >> ${path}/${output}.dnaworks.job
	echo "//" >> ${path}/${output}.dnaworks.job

	/software/DNAWorks/dnaworks ${path}/${output}.dnaworks.job

	awk '/------/ { if( p == 2 ) { exit } } { if( p == 2) { print $0 } } /The DNA sequence/ { p = 1 } /------/ { if(p == 1) { p = 2 } }'  ${path}/${output}.dnaworks.out | awk '{printf "%s", $2} END {printf "\n" }' > ${path}/${output}.seq
	echo $seq >> ${path}/${output}.seq

	/bin/rm -f ${path}/${output}.dnaworks.out ${path}/${output}.dnaworks.job
done

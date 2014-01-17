#!/bin/bash

#generate a nucmer alignment, using default parameters
nucmer --mum $REFERENCE $ASSEMBLY > nucmer.log
delta-filter -r out.delta > out.filter
show-coords -rc out.filter > out.coords

# for each chromosome
for i in X 2L 2R 3L 3R 4 XHet 2LHet 2RHet 3LHet 3RHet YHet M_iso1 U Uextra
do
printf "\n\n$i\n"

# filter out only high stringency alignments, count alignments, gaps, and total length aligned, used to generate table 2

cat out.coords | sed -e '1,5d' | awk '{if ($10>99 && $7>1000) print}' | awk -v ref="$i" '{if ($15==ref) print}' | awk '{print $16}' | sort | uniq | wc -l

cat out.coords | sed -e '1,5d' | awk '{if ($10>99 && $7>1000) print}' | awk -v ref="$i" '{if ($15==ref) print $15"\t"$1"\t"$2}' | bedtools merge | awk '{print $3-$2}' | awk '{sum+=$1} END {print sum}'

cat out.coords | sed -e '1,5d' | awk '{if ($10>99 && $7>1000) print}' | awk -v ref="$i" '{if ($15==ref) print $15"\t"$1-1"\t"$2-1}' | bedtools merge > $i.align

bedtools complement -i $i.align -g dmel-all-chromosome-r5.53.genome | awk -v ref="$i" '{if ($1==ref) print}' | wc -l

bedtools complement -i $i.align -g dmel-all-chromosome-r5.53.genome | awk -v ref="$i" '{if ($1==ref) print}'  >> gaps.bed 

done
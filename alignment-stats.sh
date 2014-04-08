#!/bin/bash

#generate a nucmer alignment, using default parameters
nucmer $REFERENCE $ASSEMBLY > nucmer.log

#filter to select only the best mapping of each contig to the reference
~/progs/MUMmer3.23/delta-filter -q out.delta > out.q.delta

#display alignments
~/progs/MUMmer3.23/show-coords -THrcl out.q.delta

# require alignments of at least 99% identity and 1000 bp
~/progs/MUMmer3.23/show-coords -THrcl out.q.delta | awk '{if ($7>99 && $5>1000) print $12"\t"$1"\t"$2"\t"$13"\t"$11}' > nucmer.bed


bedtools merge -i nucmer.bed > nucmer.merge.bed

for i in X 2L 2R 3L 3R 4 XHet 2LHet 2RHet 3LHet 3RHet YHet M U
do
    echo $i
# count the alignments
    cat nucmer.bed | awk -v i=$i '{if ($1==i) print}' | cut -f4 | sort | uniq | wc -l

# count the gaps: contiguous alignments minus 1, plus 2 for the ends of each alignment
    bedtools complement -g ../../reference/dmel-quast-euchrom-heterochrom.genome -i nucmer.merge.bed > nucmer.complement.bed
    cat nucmer.complement.bed | awk -v i=$i '{if ($1==i) print}' | wc -l

# sum the total aligned length
    cat nucmer.merge.bed | awk -v i=$i '{if ($1==i) print $3-$2}' | awk '{sum+=$1} END {print sum}'
    printf "\n\n"
done
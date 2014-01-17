#!/bin/bash

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam 2L:1-23011544 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > 2L.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam 2LHet:1-368872 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > 2LHet.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam 2R:1-21146708 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > 2R.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam 2RHet:1-3288761 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > 2RHet.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam 3L:1-24543557 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > 3L.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam 3LHet:1-2555491 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > 3LHet.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam 3R:1-27905053 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > 3R.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam M_iso1:1-17222 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > M_iso1.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam U:1-10049037 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > U.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam Uextra:1-29004656 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > Uextra.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam X:1-22422827 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > X.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam XHet:1-204112 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > XHet.mismatch

nohup java -cp sam-1.77.jar:picard-1.77.jar:. Biostar59647 dmel-all-chromosome-r5.53.fasta dmel-filter.sorted.calmd.bam YHet:1-347038 | tr '<' '\n' | sed 's/\/>//g' | sed 's/\/name//g' | sed 's/>/\t/g' | grep 'name\|mismatch' | awk '{if ($1=="name") id=$2 ; else if ($1=="M") print id"\t"$0}' | tr '"' '\t' > YHet.mismatch


#concatenate mismatches for all chromosomes
for i in *.mismatch ; do cat $i ; done >> mismatch.R.txt


#repeat above procedure for insertions and deletions
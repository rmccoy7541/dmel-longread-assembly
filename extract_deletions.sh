#!/bin/bash

# this is an incredibly inelegant script to grab the positions of deletions from the CIGAR string
# in the case of this particular BAM file, we need to go up to 127!

samtools view deletions.bam | cut -f1,6 | sed 's/S/\t/g' | sed 's/M/\t/g' | sed 's/I/\t/g' | sed 's/H/\t/g' | sed 's/D/\tD\t/g' > deletions.perl.in

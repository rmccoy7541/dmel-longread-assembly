#!/bin/bash

# this script extracts the positions of deletions in a read by reading from the CIGAR string
# assumes that deletions.bam includes only those reads that contain a deletion

samtools view deletions.bam | cut -f1,6 | sed 's/S/\t/g' | sed 's/M/\t/g' | sed 's/I/\t/g' | sed 's/H/\t/g' | sed 's/D/\tD\t/g' > deletions.perl.in

cat deletions.perl.in | perl generate_deletion_profile.pl > deletions.positions

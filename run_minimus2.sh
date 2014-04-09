#!/bin/bash

# convert celera assembly to Amos format
toAmos -s celera_output.fa -o minimus_input.afg

# run minimus2 on celera assembly
minimus2 minimus_input -D REFCOUNT=0 -D MINID=99.9 -D OVERLAP=800 -D MAXTRIM=1000 -D WIGGLE=15 -D CONSERR=0.01

# concatenate supercontigs and singletons (i.e. celera contigs)
cat minimus_input.fasta minimus_input.singletons.seq > minimus_output.fa

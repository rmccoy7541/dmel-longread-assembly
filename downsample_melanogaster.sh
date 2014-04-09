#!/bin/bash

# Start processing options at index 1.
OPTIND=1
# OPTERR=1
while getopts ":c:d:o:f:" VALUE "$@" ; do
    if [ "$VALUE" = "c" ] ; then
		COVERAGE=$OPTARG
    fi
    if [ "$VALUE" = "d" ] ; then
		DATA_FILE=$OPTARG
    fi
    if [ "$VALUE" = "o" ] ; then
		OUT_DIR=$OPTARG
    fi
    if [ "$VALUE" = "f" ] ; then
		OUT_FILE=$OPTARG
    fi
    # The getopt routine returns a colon when it encounters
    # a flag that should have an argument but doesn't.  It
    # returns the errant flag in the OPTARG variable.
    if [ "$VALUE" = ":" ] ; then
        echo "Flag -$OPTARG requires an argument."
        echo "Usage: $0 [-x coverage] [-d data directory] \
			[-o output file]"
        exit 1
    fi
    # The getopt routine returns a question mark when it
    # encounters an unknown flag.  It returns the unknown
    # flag in the OPTARG variable.
    if [ "$VALUE" = "?" ] ; then
        echo "Unknown flag -$OPTARG detected."
        echo "Usage: $0 [-x coverage] [-d data directory] \
			[-o output file]"
        exit 1
    fi
done

# Get total reads in raw, unsampled fq
READLINES=`wc -l $DATA_FILE | cut -d " " -f 1`

# 4 lines per read in fq file

READS=`python -c "print(int('$READLINES')/4)"`

# These reads = 34x coverage so divide by 34 then multiply by coverage

# Downsample to X coverage
DSREADS=`python -c "print(int(round(((int('$READS')/34) * float('$COVERAGE')))))"`

seqtk sample $DATA_FILE $DSREADS > $OUT_DIR/$OUT_FILE.fq

# Now make frg for wgs assembly.
cd $OUT_DIR
/home/rwtaylor/bin/wgs-8.1/Linux-amd64/bin/fastqToCA \
-libraryname $OUT_FILE \
-technology 'moleculo' \
-type 'sanger' \
-reads $OUT_DIR/$OUT_FILE.fq \
> $OUT_FILE.frg


# Unforce BOHunitigger
sed -i 's/forceBOGunitigger=1/forceBOGunitigger=0/' $OUT_FILE.frg


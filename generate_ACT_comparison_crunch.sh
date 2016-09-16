#!/bin/bash

# follow procedures in supplementary data of this paper: http://www.nature.com/nprot/journal/v7/n7/full/nprot.2012.068.html

if [ $# -ne 4 ]; then 
    echo ""
    echo "Usage: $0 [blastn/tblastx] [input_first_genome.fna] [input_second_genome.fna] [output_filename]"
    echo ""
    exit 1
fi



#formatdb="/vol/biotools/bin/formatdb"
#blastall="/vol/biotools/bin/blastall"

# run formatdb
formatdb -p F -i $2 

# do blastall using 30 cores
blastall -p $1 -m 8 -e 1e-5 -d $2 -i $3 -o $4 -a 30

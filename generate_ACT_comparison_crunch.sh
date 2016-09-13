#!/bin/bash

# follow procedures in supplementary data of this paper: http://www.nature.com/nprot/journal/v7/n7/full/nprot.2012.068.html

if [ $# -ne 3 ]; then 
    echo ""
    echo "Usage: $0 [input_first_genome.fna] [input_second_genome.fnaa] [output_filename]"
    echo ""
    exit 1
fi



#formatdb="/vol/biotools/bin/formatdb"
#blastall="/vol/biotools/bin/blastall"

# run formatdb
formatdb -p F -i $1 

# do blastall using 30 cores
blastall -p blastn -m 8 -e 1e-5 -d $1 -i $2 -o $3 -a 30

#!/bin/bash

if [ $# != 2 ]; then
    echo ""
    echo "Usage: $0 [input_bam_file] [input_interested_contigs.fasta]"
    echo ""
    exit 1
fi


BAMFILE=$1
CONTIGS=$2
OUTPUT=${CONTIGS%.*}

#echo "=> slice bam file by input contigs in fasta format, generate sliced bam/readname files"
grep -e ">" $CONTIGS | awk 'sub(/^>/, "")' > $OUTPUT.contig
samtools view -bh $BAMFILE `cat $OUTPUT.contig | tr "\n" " "` > $OUTPUT.bam
samtools view $OUTPUT.bam | awk '{print $1}' # > $OUTPUT.readname
#echo "=> done!"

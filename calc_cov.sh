#!/bin/bash


# ---------
# Software and function
# ---------
java=/home/hou/Software/java/java-8-openjdk-amd64/bin/java
bwa="/home/hou/Software/Module_Alignment/bwa/bwa"
samtools="/home/hou/Software/Module_NGS_tools/samtools/samtools"


function bwa_align_paired_reads {
    # test input  
    if [ "$#" -ne 4 ]; then
        echo "Error in bwa_align_paired_reads, input arguments should be: [reference.fasta, fastq1, fastq2, output_basen
ame]!"
        exit 1
    fi

    # get file basename, also remove the given suffix, here fasta
    filename=$(basename "$1")
    extension="${filename##*.}"
    filestem="${filename%.*}"

    #Index reference
    if [[ ! -e ${1}.bwt ]]; then
        echo -e "run CMD: '$bwa index $1'"
        $bwa index $1
    fi

    #Sort reference
    if [[ ! -e ${1}.fai ]]; then
        echo -e "run CMD: '$samtools faidx $1'"
        $samtools faidx $1
    fi

    echo -e "run CMD: $bwa mem -t 16 $1 $2 $3 > $4_paired.sam"
    $bwa mem -t 30 $1 $2 $3 > $4_paired.sam
    $samtools view -bS $4_paired.sam | $samtools sort -@ 30 - $4_sorted
    $samtools index $4_sorted.bam
    echo -e "done\n"
    }


function get_coverage_from_bam {
    # test input  
    if [ "$#" -ne 2 ]; then
        echo "Usage: $0: [input_bam_file] [output_prefix]"
        exit 1
    fi

    # do samtools depth
    $samtools depth $1 > $2_depth.txt

    # do calcualte coverage
    python ./calc_cov_in_bam_depth.py -i $2_depth.txt -o $2_coverage.txt
       
}



# ---------
# arguments
# ---------

if [ $# != 4 ]; then
    echo -e "\nUsage: $0 [input_reference] [input_fwd_fastq] [input_rev_fastq] [prefix]\n"
    exit 1
fi



# ---------
# alignment
# ---------

bwa_align_paired_reads $1 $2 $3 $4

# --------
# bam2coverage
# --------
get_coverage_from_bam $4_sorted.bam $4_sorted


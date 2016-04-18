#!/usr/bin/env python


# Copyright (C) 2016  Shengwei Hou : housw2010 'at' gmail 'dot' com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import sys
import subprocess
import os
import argparse
from utils.Fastq import parse_fastq


def write_fastq_records(input_fwd_fastq, input_rev_fastq, fastq_header_set, prefix):
    """
    :param input_fwd_fastq: input forward reads in fastq format
    :param input_rev_fastq: input reverse reads in fastq format
    :param fastq_header_set: input a set of fastq headers
    :param prefix: output prefix of extracted fastq
    :return: Nonething

    @FCC57ARACXX:1:1103:6389:21631#/2
    AAACCACCACTTCCTAATTCTAATTGTAAATGTGAAAATTGTCTCACAATTTCATTACGAATGGTATTTACTGCTGCTGCATCTCGAAGTGTTGGTTTAA
    +
    a_aecce`egggghddffhfdffidhcbhfhhgb`gfhhhfggdgiihiifggghihihgfheeWaeghhhhhhffffhi_geggb_aaZ_bc^a_aabb
    """
    for fastq_file in (input_fwd_fastq, input_rev_fastq):
        basename = os.path.basename(fastq_file)
        output_file = prefix + basename

        with open(output_file, "w") as oh:
            fastq_records = parse_fastq(fastq_file)
            for record in fastq_records:
                name = record.header.rstrip("/1").rstrip("/2")
                #print name
                if name in fastq_header_set:
                    oh.write(str(record))


def get_fastq_headers(input_bam, input_contigs):
    """
    :param input_bam:     input bam file contains read alignment to contigs
    :param input_contigs: a subset of interested contigs
    :param prefix: output prefix
    :return:  a set of fastq headers
    """
    fastq_headers = []

    p = subprocess.Popen(["bash", "slice_bam_by_contigs.sh", input_bam, input_contigs],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    for header in stdout.split("\n"):
        fastq_headers.append(header.strip())

    return set(fastq_headers)


def main():

    # parse arguments
    parser = argparse.ArgumentParser(description="replace spaces in fasta header by underscores")
    parser.add_argument("input_bam", help="input bam file contains read alignments to contigs")
    parser.add_argument("interest_contigs", help="a subset of interested contigs in fasta format to fish the raw reads")
    parser.add_argument("fwd_fastq", help="the forward fastq reads")
    parser.add_argument("rev_fastq", help="the reverse fastq reads")
    parser.add_argument("-p", "--prefix", required=False, default="extracted_", help="the output prefix of extracted fastq files")

    args = parser.parse_args()

    # get headers
    fastq_header_set = get_fastq_headers(args.input_bam, args.interest_contigs)

    # write fastq for given headers
    write_fastq_records(args.fwd_fastq, args.rev_fastq, fastq_header_set, args.prefix)







if __name__ == "__main__":
    main()
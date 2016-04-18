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
import os
import argparse
from utils.Fasta import parse_fasta
from utils.SAM import SamParser
from utils.Fastq import parse_fastq


def get_interest_contig_headers(input_contigs):
    """
    :param input_contigs: input contigs file in fasta format
    :return: a set of contig headers
    """
    headers = []

    contigs = parse_fasta(input_contigs)
    for contig in contigs:
        headers.append(contig.header)

    return set(headers)


def get_fastq_headers(input_sam, contig_header_set):
    """
    :param input_sam: input sam file contains read alignments to contigs
    :param contig_header_set: interested contigs to subset
    :return:  a set of fastq headers, which aligned to interested contigs

    FCC57ARACXX:1:1310:12325:8819#	99	Contig_31	6976	42	100M	=	7044	168	CATTGTTGCCAAGAGGTCCTAGTCCCTGAC
    CGCGGATCCTCACTGTTGCCGTCGAACTTGGATTACCTGTGGTGGTGATTTGTAAACCAGGCACTCTTCC	_bbeeeeegegfgiii_fghffdedgddfhfiifhaf`cgifhhiiifhihggcacdccd`b
    bbcccbb]`^[^a\W___ccccYab[Wa^aaabX_b]_	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:100	YS:i:0	YT:Z:CP

    """
    headers = []

    sam_records = SamParser(input_sam)
    for sam in sam_records:
        if sam.rname in contig_header_set:
            headers.append(sam.rname)

    return set(headers)


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
                if record.header.rstrip("/1").rstrip("/2") in fastq_header_set:
                    oh.write(str(record))




def main():

    # parse arguments
    parser = argparse.ArgumentParser(description="replace spaces in fasta header by underscores")
    parser.add_argument("input_sam", help="input sam file contains read alignments to contigs")
    parser.add_argument("interest_contigs", help="a subset of interested contigs in fasta format to fish the raw reads")
    parser.add_argument("fwd_fastq", help="the forward fastq reads")
    parser.add_argument("rev_fastq", help="the reverse fastq reads")
    parser.add_argument("-p", "--prefix", required=False, default="extracted_", help="the output prefix of extracted fastq files")

    args = parser.parse_args()
    input_sam = args.input_sam
    input_contigs = args.interest_contigs
    fwd_fastq = args.fwd_fastq
    rev_fastq = args.rev_fastq
    prefix = args.prefix

    print input_sam, input_contigs, fwd_fastq, rev_fastq, prefix

    contig_header_set = get_interest_contig_headers(input_contigs)
    fastq_header_set = get_fastq_headers(input_sam, contig_header_set)
    write_fastq_records(fwd_fastq, rev_fastq, fastq_header_set, prefix)







if __name__ == "__main__":
    main()
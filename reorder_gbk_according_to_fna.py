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
from Bio import SeqIO


def parse_genbank_block(gb_file):
    """ parse a block of genbank record, which follows the following pattern:
        
        1) starts with a "LOCUS"" keyword, followed by spaces, then seqname
        2) ends with "//" sign 
        3) contents in between will not be parsed

        ^LOCUS        seqname length bp DNA linear date
        ...
        ...
        ^//

        this function will yield each genbank record as a generator

    """
    seqname = None
    block = [] # use list instead of string to improve performance

    with open(gb_file, "r") as ih:
        for line in ih:
            # a new genbank record 
            if line.startswith("LOCUS"):
                # the previous block wasn't end with a "//", yield with a Warning
                if block:
                    print "WARNING: The previous block wasn't ended with a '//'!!!\n"
                    yield (seqname, block)
                line_list = line.strip().split(" ")
                line_list = [item for item in line_list if item != ""] # remove "" items
                seqname = line_list[1]
                block = [line]
            
            # the end of a genbank record
            elif line.startswith("//"):
                block += line
                yield (seqname, "".join(block))
                block = []
                seqname = None

            # lines within a genbank record
            else:
                block.append(line)



def main():
    parser = argparse.ArgumentParser(description="order genbank file according to fna file")
    parser.add_argument("-f", "--fna", required=True, help="input fna file")
    parser.add_argument("-g", "--genbank", required=True, help="input genbank file")
    parser.add_argument("-p", "--prefix", required=False, help="output prefix")
    args = parser.parse_args()

    # test prefix 
    if not args.prefix:
        basename = os.path.basename(args.genbank)
        prefix = os.path.splitext(basename)[0]
    else:
        prefix = args.prefix
    
    # read in fna file, save record names
    contig_names = []
    for rec in SeqIO.parse(args.fna, "fasta"):
        contig_names.append(rec.name)
    
    # read in genbank file, get seqname:genbank dict
    seqname2genbank = {} # {seqname:genbank}
    for seqname, rec in parse_genbank_block(args.genbank):
        print seqname + " was parsed from genbank file!"
        seqname2genbank[seqname] = rec

    # write out genbank record in order given by contig_names
    with open(prefix+"_reordered.gbk", "w") as oh:
        for seqname in contig_names:
            assert seqname in seqname2genbank, "seqname: %s in fna wasn't found in genbank file!!!"%seqname
            genbank = seqname2genbank[seqname]
            oh.write(genbank)
 
if __name__ == '__main__':
    main()
    
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


import gzip


class Fastq(object):

    def __init__(self, header, seq, qual):
        self.header = header
        self.seq = seq.upper()
        self.qual = qual

    def print_pretty_fasta(self, width=100):
        """
        :param width: maximum 'width' per line for sequences
        :return:      return pretty string representation of Fasta record
        """
        ret = ">" + self.header + "\n"
        for i, char in enumerate(self.seq):
            if i > 0 and i % width == 0:
                ret += "\n"
            ret += char
        return ret + "\n"

    def __repr__(self):
        return 'Fastq(header=%s, seq=%s, qual=%s)' % (self.header, self.seq, self.qual)

    def __str__(self):
        return "@" + self.header + "\n" \
               + self.seq +"\n" \
               + "+\n" \
               + self.qual + "\n"


def parse_fastq(fastq_file):
    """
    :param fastq_file: input fastq file
    :return:           yield Fastq record as a generator
    """

    # handle gzip compressed file
    if fastq_file.endswith(".gz"):
        ih = gzip.open(fastq_file, "r")
    else:
        ih = open(fastq_file, "r")

    # parse fastq records
    header = ih.readline().strip()[1:]
    while header:
        header = header.strip()[1:]
        seq = ih.readline().strip()
        ih.readline()
        qual = ih.readline().strip()
        yield Fastq(header, seq, qual)
        header = ih.readline()

    ih.close()
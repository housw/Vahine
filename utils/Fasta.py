#!/usr/bin/env python


# <Fasta.py, utilities for fasta format manipulations>
# Copyright (C) <2016>  <Shengwei Hou> <housw2010'at'gmail'dot'com>
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

rev_dict = {"A": "T",
            "T": "A",
            "C": "G",
            "G": "C"}


class Fasta(object):

    def __init__(self, header, seq):
        self.header = header
        self.seq = seq.upper()

    def get_gc_content(self):
        ret = 0
        if len(self.seq) == 0:
            return ret
        for c in self.seq:
            if c == "G" or c == "C":
                ret += 1
        return float(ret)/len(self.seq)

    def __str__(self):
        return ">"+self.header+"\n"+self.seq+"\n"


def parse_fasta(fasta_file):
    """
    :param fasta_file: input fasta file
    :return:           yield Fasta record as a generator
    """
    header = ""
    seq = []
    with open(fasta_file, "r") as ih:
        for line in ih:
            if line.startswith(">"):
                if header:
                    yield Fasta(header, "".join(seq))
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        yield Fasta(header, "".join(seq))










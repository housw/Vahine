#!/usr/bin/env python


# <Blast.py, utilities for Blast format manipulations>
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


class BlastRecord(object):
    """ This class used to represent tabular blast record in *.tab file, fields are

        0     1     2      3         4      5     6     7      8      9      10     11
    QID RefID Ident AlignLen MisMatch Gaps QStart QEnd RefStart RefEnd Evalue BitScore

    """

    __slots__=["QID", "RefID", "Ident", "AlignLen", "MisMatch", "Gaps",
               "QStart", "QEnd", "RefStart", "RefEnd", "Evalue", "BitScore"]

    def __init__(self,
                 QID=None, RefID=None,
                 Ident=None, AlignLen=None,
                 MisMatch=None, Gaps=None,
                 QStart=None, QEnd=None,
                 RefStart=None, RefEnd=None,
                 Evalue=None, BitScore=None
                 ):
        self.QID=QID
        self.RefID=RefID
        self.Ident=Ident
        self.AlignLen=AlignLen
        self.MisMatch=MisMatch
        self.Gaps=Gaps
        self.QStart=QStart
        self.QEnd=QEnd
        self.RefStart=RefStart
        self.RefEnd=RefEnd
        self.Evalue=Evalue
        self.BitScore=BitScore

    def __str__(self):
        """ string representation of blast record, the same with initial blast
            tabular output
        """
        _str=""
        _str += self.QID+"\t"
        _str += self.RefID+"\t"
        _str += self.Ident+"\t"
        _str += str(self.AlignLen)+"\t"
        _str += str(self.MisMatch)+"\t"
        _str += str(self.Gaps)+"\t"
        _str += str(self.QStart)+"\t"
        _str += str(self.QEnd)+"\t"
        _str += str(self.RefStart)+"\t"
        _str += str(self.RefEnd)+"\t"
        _str += str(self.Evalue)+"\t"
        _str += str(self.BitScore)
        return _str


class BlastRecordParser(object):
    """ This class used to parse tabular blast record in *.tab file, yield
        BlastRecord instance one by one
    """

    def __init__(self, fh):
        self.handle = fh

    def _yield_records(self):
        """ This function used to yield tabular blast records
        """
        _file_handle = None
        if not hasattr(self.handle, "read"):
            try:
                _file_handle = open(self.handle, "r")
            except Exception as e:
                print e
                print "can not open file for read !"
        else:
            _file_handle = self.handle

        for i, line in enumerate(_file_handle):
            if line.startswith("\n") or line.startswith("#"):
                continue
            line = line.strip().split("\t")
            assert len(line) == 12, "Bad blast record was found in %d-th line"%i+1
            yield BlastRecord(*line)

    def __iter__(self):
        return self._yield_records()
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


class CentrifugeRecord(object):
    """ This class used to represent centrifuge records in tabular result

        0      1       2     3         4           5        6
    readID uniqueID  taxID score secBestScore hitLength numMatches

    """
    __slots__=["readID", "uniqueID", "taxID", "score", "secBestScore",
               "hitLength", "numMatches"]

    record_dict = {} # {readID:CentrifugeRecord}

    def __init__(self, readID, uniqueID, taxID, score, secBestScore, hitLength, numMatches):
        self.readID = readID
        self.uniqueID = uniqueID
        self.taxID = taxID
        self.score = score
        self.secBestScore = secBestScore
        self.hitLength = hitLength
        self.numMatches = numMatches

        if not self.readID in CentrifugeRecord.record_dict:
            CentrifugeRecord.record_dict[self.readID] = [self]
        else:
            CentrifugeRecord.record_dict[self.readID].append(self)

    def __str__(self):
        return self.readID +"\t"+self.uniqueID+"\t"+\
               str(self.taxID)+"\t"+str(self.score)+"\t"+str(self.secBestScore)+\
               "\t"+str(self.hitLength)+"\t"+str(self.numMatches)


class CentrifugeRecordParser(object):
    """ parse centrifuge results as a generator
    """
    def __init__(self, fh):
        """This parser is used to parse centrifuge result
        """
        self.handle = fh

    def _yield_records(self):
        """ This function used to yield Centrifuge record one by one
        """
        _file_handle = None
        # a opened file
        if hasattr(self.handle, 'read'):
            _file_handle = self.handle
        else:
            try:
                _file_handle = open(self.handle, "r")
            except Exception as e:
                print "cann't open file for read: %s !"%e

        for i, line in enumerate(_file_handle):
            if line.startswith("readID") or \
               line.startswith("\n") or \
               line.startswith(" ") or \
               line.startswith("\t"):
                continue
            line = line.strip().split("\t")
            #assert len(line) == 7, "The lenght of %d-th line is not 7 !!!"
            if len(line) != 7:
                print "Warning: the fields number of line %d is not 7, will be passed!!"%i
                print line
                continue
            yield CentrifugeRecord(*line)

    def __iter__(self):
        return self._yield_records()
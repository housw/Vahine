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


class SamRecord(object):
    """ class to represent the SAM record in segemehl sam file
    """

    __slots__ = ["qname","shortname", "oriStrand","weight", "mapStrand",
                 "rname", "startPos", "mapq", "cigar",
                 "rnext", "pnext", "tlen", "seq","qual", "option"]

    def __init__(self, qname, shortname, oriStrand, weight, mapStrand, rname,
                 startPos, mapq, cigar, rnext, pnext, tlen,
                 seq, qual, option):

        self.qname = qname
        self.shortname = shortname
        self.oriStrand = oriStrand
        self.weight = weight
        self.mapStrand = mapStrand
        self.rname = rname
        self.startPos = startPos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        self.option = option
        #self._update_weight()

    def _update_weight(self):
        """ divide the weight according to # of hit positions
        """
        _repeat_time = int(self.option['NH'])
        self.weight /= float(_repeat_time)

    def __str__(self):
        return self.qname + " mapped to "+self.mapStrand + \
        " strand, locates in <" + self.rname + "> starts at <" + \
        str(self.startPos)+">\n"



class SamParser(object):
    """ class to parse the segemehl SAM record from a open file handle
    """

    def __init__(self, handle_or_fileStr):

        self.handle = handle_or_fileStr

    def _lines(self):

        # judge file opened or not
        if not hasattr(self.handle, "read"):
            handle = open(self.handle, "r")
        else:
            handle = self.handle

        while True:
            line = handle.readline()
            if not line:
                handle.close()
                break
            else:
                if line.startswith("@"):
                    continue
                else:
                    line = line.strip()
                    yield line

    def _parse_SAM(self, inline):
        for line in inline:
            sam_line = line.strip().split("\t")
            fullname = sam_line[0].replace("_", " ")
            qname = " ".join(item for item in fullname.split(" ")[0:2])
            shortname = fullname.split(" ")[0]
            if " " in fullname:
                oriStrand = int(fullname.split(" ")[1][0])
            else:
                oriStrand = None
            if "weight|" in fullname:
                weight = float(fullname.strip().split("|")[1])
            else:
                weight = 1

            # the strand info that reads mapped in reference genome
            if int(sam_line[1]) == 16 or int(sam_line[1]) == 272:
                mapStrand = "-"
            elif int(sam_line[1]) == 0 or int(sam_line[1]) == 256:
                mapStrand = "+"
            else:
                mapStrand = None

            # cigar info and use it to get map info
            cigar = sam_line[5]
            # use left most coordinate system
            startPos = int(sam_line[3])
            rname = sam_line[2]
            mapq = int(sam_line[4])
            rnext = sam_line[6]
            pnext = int(sam_line[7])
            tlen = int(sam_line[8])
            seq = sam_line[9]
            qual = sam_line[10]

            # put all additional fields into option
            option = {}

            for item in sam_line[11:]:
                parts = item.strip().split("\t")
                for part in parts:
                        new_part = part.strip().split(":")
                        option.update({new_part[0]: new_part[2]})

            sam_record = SamRecord(qname, shortname, oriStrand, weight, mapStrand,
                                   rname, startPos, mapq, cigar,
                                   rnext, pnext, tlen, seq, qual, option)

            yield sam_record

    def __iter__(self):
        return self._parse_SAM(self._lines())
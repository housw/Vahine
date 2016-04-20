#!/bin/bash

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


if [ $# != 2 ]; then
    echo ""
    echo "Usage: $0 [output_dir] [out_filename]"
    echo ""
    exit 1
fi

ALL=ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz
OUT_DIR=$1
FILENAME=$2


# ----------
# Download all.fna.tar.gz, and extract all folders
# ----------
if [ ! -f "$OUT_DIR" ]; then
    mkdir $OUT_DIR
fi

echo "Start to download $ALL"
wget $ALL -O $OUT_DIR/all.fna.tar.gz 
tar zxf $OUT_DIR/all.fna.tar.gz -C $OUT_DIR
rm $OUT_DIR/all.fna.tar.gz


# -----------
# write out all fna to out_filename
# -----------
if [ -f "$FILENAME" ]; then
    echo "File $FILE exists, remove it!"
    rm $FILENAME
fi

for f in `ls $OUT_DIR`; do
    #echo "$f"

    for nc in `ls $OUT_DIR/$f`; do
        if [[ ! -z $(echo ${nc} | grep fna) ]]; then
            echo "writing $OUT_DIR/$f/$nc"
            cat $OUT_DIR/$f/$nc >> $FILENAME
        fi 
    done  
done

echo "Done!"

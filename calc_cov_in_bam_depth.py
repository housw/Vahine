#!/usr/bin/python

import sys, os 
import argparse


def main():
    parser = argparse.ArgumentParser(description="Use this script to calculate coverage of each contig in a samtools depth file.")
    parser.add_argument('-i', '--input_depth', required=True, action="store", dest="depth_file")
    parser.add_argument('-o', '--output_file', required=True, action="store", dest="output_file")
    args = parser.parse_args()

    coverage = {}
    order = []
    length = {}
    effective_length = {}

    # open depth file 
    with open(args.depth_file, "r") as ih:
        for line in ih:
            splitline = line.strip().split("\t")
            contig, pos, cov = splitline[0], int(splitline[1]), int(splitline[2])
            # if not start of contig
            if contig in coverage:
                coverage[contig] += cov
                length[contig] = pos
                if cov >= 1:
                    effective_length[contig] += 1
            # new contig starts
            else:
                coverage[contig] = cov
                length[contig] = pos
                if cov >= 1:
                    effective_length[contig] = 1
                else:
                    effective_length[contig] = 0
                # remember input contig order
                order.append(contig) 

    # write out depth file 
    with open(args.output_file, "w") as oh:
        oh.write("#Name\tSumCov\tAvgCov\tRefLen\tEffCov\tEffLen\tCovPct\n")
        for contig in order:
            SumCov = coverage[contig]
            RefLen = length[contig]
            EffLen = effective_length[contig]
            AvgCov = float(SumCov)/RefLen
            EffCov = float(SumCov)/EffLen
            CovPct = 100*float(EffLen)/RefLen
            oh.write(contig+"\t"+str(SumCov)+"\t%.2f\t"%AvgCov+str(RefLen)+"\t%.2f\t"%EffCov+str(EffLen)+"\t%.2f\n"%CovPct)
            


if __name__ == "__main__":
    main()

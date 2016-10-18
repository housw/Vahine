#!/usr/bin/env python

import os, sys
import argparse
from ftplib import FTP
import ftputil
import subprocess


def main():

    description="""
    Select and filter your interest genomes from NCBI genome: https://www.ncbi.nlm.nih.gov/genome/browse/,
    then download the assembly reports of selected records in tab-deliminated txt file.
    Download all the genome sequences in the downloaded assembly reports.
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', "--input_file", type=str, required=True, help=description)
    parser.add_argument("-o", "--outdir", type=str, help="output folder, default=genomes", default="genomes")
    args = parser.parse_args()

    # display help
    if len(sys.argv) == 1:
        print parser.format_help()
        sys.exit(1)


# ----------
#  parse assembly report to get assembly ID and ftp path
# ----------

    assembly2ftp = {} # {assembly: ftp}
    with open(args.input_file, "r") as ih:
        header = ih.readline().strip().split('\t')
        assemby_idx = header.index('Assembly')
        ftp_idx = header.index('GenBank FTP')
        for line in ih:
            line = line.strip().split("\t")
            assembly = line[assemby_idx].strip()
            ftp = line[ftp_idx].strip()
            assert assembly not in assembly2ftp, "Repeating genbank assembly ID %s has been found!"%assembly
            assembly2ftp.update({assembly:ftp})



# --------------------
# download assemblies
# --------------------
    # login NCBI ftp site
    host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'housw2010@gmail.com')

    # download all the files in the directory
    for assembly, ftp_path in assembly2ftp.items():

        # create a local folder
        local_folder = os.path.join(args.outdir, assembly)
        subprocess.check_call(['mkdir', '-p', local_folder])

        # move to the ftp_path folder on host
        server_folder = ftp_path.split("ftp.ncbi.nlm.nih.gov")[-1]
        host.chdir(server_folder)
        file_list = host.listdir(host.curdir)
        for remote_file in file_list:
            if host.path.isfile(remote_file):
                print "Start to download %s ..."%remote_file
                local_file = os.path.join(local_folder, remote_file)
                # download arguments: remote source filename, local target filename, callback
                host.download(remote_file, local_file)


if __name__ == "__main__":
    main()

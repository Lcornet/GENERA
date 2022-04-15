#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""chaneg IDs in spades output"""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2021, University of Liège"
__version__ = "1.0.0"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"

import glob
import re

def main():

    #open OGs to define core gene
    fna_list = glob.glob("*.fna")
    for fna in fna_list:
        #define
        file = fna.replace(".fna", "")
        name = str(file) + '-spades.fna'
        out_file = open((name), "w")
        fnafile = open(fna)
        for line in fnafile:
            file_record = line.replace("\n", "")
            if '>' in line:
                defline = re.sub(r'\_length_\d+','', file_record)
                defline = re.sub(r'\_cov_\d+','', defline)
                defline = re.sub(r'\.\d+','', defline)
                defline = re.sub(r'\_ID_\d+','', defline)
                out_file.write(str(defline) + "\n")
            else:
                out_file.write(str(file_record) + "\n")

if __name__ == '__main__':
    main()
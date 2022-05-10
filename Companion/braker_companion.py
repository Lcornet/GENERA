#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Kraken-ids-creation.py: Creation of ids to add genome in kraken db library."""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2020, University of Liège"
__version__ = "1.0.0"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"

import click
import glob

def main():

    out = open("bam_list.txt", "w")
    bam_list = glob.glob("*.bam")
    length = len(bam_list)
    num = 0
    for bam in bam_list:
        num += 1
        if (num < length):
            out.write(str(bam) + ",")
        else:
            out.write(str(bam) + "\n")

if __name__ == '__main__':
    main()
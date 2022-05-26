#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""change-ext.py: Change ext of fasta/ali files."""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2020, University of Liège"
__version__ = "1.0.1"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"

import click
import glob

@click.command()
###ARGUMENT#####
#taxid
@click.option('--currentext', default='fasta', help='mode')
@click.option('--newext', default='faa', help='mode')

def main(currentext, newext):

    if (currentext == 'fasta'):
        fasta_list = glob.glob("*.fasta")

        for fasta in fasta_list:
            basename = fasta.replace(".fasta", "")
            mod_name = str(basename) + "." + newext
            mod_file = open(mod_name, "w")
            fastafile = open(fasta)
            for line in fastafile:
                record = line.replace("\n", "")
                mod_file.write(str(record) + "\n")
    
    if (currentext == 'fa'):
        fasta_list = glob.glob("*.fa")

        for fasta in fasta_list:
            basename = fasta.replace(".fa", "")
            mod_name = str(basename) + "." + newext
            mod_file = open(mod_name, "w")
            fastafile = open(fasta)
            for line in fastafile:
                record = line.replace("\n", "")
                mod_file.write(str(record) + "\n") 

    if (currentext == 'ali'):
        fasta_list = glob.glob("*.ali")

        for fasta in fasta_list:
            basename = fasta.replace(".ali", "")
            mod_name = str(basename) + "." + newext
            mod_file = open(mod_name, "w")
            fastafile = open(fasta)
            for line in fastafile:
                record = line.replace("\n", "")
                mod_file.write(str(record) + "\n")  

    if (currentext == 'faa'):
        fasta_list = glob.glob("*.faa")

        for fasta in fasta_list:
            basename = fasta.replace(".faa", "")
            mod_name = str(basename) + "." + newext
            mod_file = open(mod_name, "w")
            fastafile = open(fasta)
            for line in fastafile:
                record = line.replace("\n", "")
                mod_file.write(str(record) + "\n") 

    if (currentext == 'fna'):
        fasta_list = glob.glob("*.fna")

        for fasta in fasta_list:
            basename = fasta.replace(".fna", "")
            mod_name = str(basename) + "." + newext
            mod_file = open(mod_name, "w")
            fastafile = open(fasta)
            for line in fastafile:
                record = line.replace("\n", "")
                mod_file.write(str(record) + "\n") 
		

if __name__ == '__main__':
    main()
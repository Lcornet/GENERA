#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""anvio_pan-to-OGs.py - extract Orthologous groups from pangennomics of anvi'o"""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2020, University of Liège"
__version__ = "1.0.0"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"

import click

@click.command()
###ARGUMENT#####
#main file - depending of the mode can be checkm result or rnammer output
@click.argument('main_file', type=click.Path(exists=True,readable=True))

def main(main_file):
    #Read main tab file
    infile_main = open(main_file)

    OGs_of = {}

    for line in infile_main:
        main_record = line.replace("\n", "")

        if ('>' in main_record):
            split_list = main_record.split("|")
            name = split_list[1]
            name = name.replace("gene_cluster:", "")
            OGs_of[name]=1
    
    #Loop in name and print OG
    for name in OGs_of:
        concat = str(name) + '_OG.fasta'
        out_file = open(concat, "w")

        infile_loop = open(main_file)
        tag = 0
        for line in infile_loop:
            main_record = line.replace("\n", "")
            if ('>' in main_record):
                split_list = main_record.split("|")
                ID = split_list[1]
                ID = ID.replace("gene_cluster:", "")        
                if (ID == name):
                    tag = 1
                else:
                    tag = 0
            
            if (tag == 1):
                out_file.write(str(main_record) + "\n")



if __name__ == '__main__':
    main()
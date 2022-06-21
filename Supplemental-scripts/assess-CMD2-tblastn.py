#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""."""
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
@click.argument('main_file', type=click.Path(exists=True,readable=True))
#taxid
@click.option('--tblastn', default='fasta', help='mode')

def main(main_file, tblastn):

    #get length of prot
    fasta_file = open(main_file)
    pastID = 'NA'
    seq = 'NA'
    length_of = {}
    for line in fasta_file:
        record = line.replace("\n", "")
        if '>' in line:
            defline = record.replace(">", "")
            split_list = defline.split(" ")
            ID = split_list[1]
            if ((ID != pastID) and (pastID != 'NA')):
                length = len(seq)
                length_of[pastID] = length
        else:
            if (seq == 'NA'):
                seq = str(record)
            else:
                seq = str(seq) + str(record) 
    length = len(seq)
    length_of[pastID] = length

    #parse tblastn
    tblastn_file = open(tblastn)
    for line in tblastn_file:
        record = line.replace("\n", "")
        split_list = defline.split("\t")
        ID = split_list[0]
        hit = split_list[1]
        ident = split_list[2]
        hitL = split_list[3]

        #Compute
        length = length_of[ID]

        if (float(ident) > 95):
            print(str(ID) + "\t" + str(length) + "\t" + str(hitL) + "\t" + str(ident) + "\t" + str(hit) + "\n")




#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Companion"""
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
@click.argument('tax_file', type=click.Path(exists=True,readable=True))
@click.argument('list_file', type=click.Path(exists=True,readable=True))

def main(tax_file, list_file):

    #Load tax file
    last_ofs = {}
    taxfile = open(tax_file)
    for line in taxfile:
        record = line.replace("\n", "")
        split_list = record.split("\t")
        genome = split_list[0]
        lineage = split_list[3]
        split_list = lineage.split("; ")
        last = split_list.pop()
        last = last.replace(" ", "-")
        last_ofs[genome]= last

    #Print according to list
    listfile = open(list_file)
    out_file = open('idm.txt', "w")
    for line in listfile:
        record = line.replace("\n", "")
        split_list = record.split("\t")
        genome = split_list[0]
        last = 'NA'
        if (genome in last_ofs):
            last = last_ofs[genome]
        else:
            last = 'NOT'
        out_file.write(str(last) + "\t" + str(genome) + "\n")

if __name__ == '__main__':
    main()




#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""deflines-rename-ali.py"""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2021, Sciensano"
__version__ = "1.0.0"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"

import click
import glob

@click.command()
###ARGUMENT#####
#main file - depending of the mode can be checkm result or rnammer output
@click.argument('main_file', type=click.Path(exists=True,readable=True))


def main(main_file):

    #Open main file, should be an IDM (long - short)
    idm_ofs = {}
    idm_file = open(main_file)
    for line in idm_file:
        record = line.replace("\n", "")
        split_list = record.split("\t")
        long = split_list[0]
        short = split_list[1]
        idm_ofs[long]=short
    
    #Make a list of ali-fasta file
    ali_list = glob.glob("*.ali")

    for ali in ali_list:
        basename = ali.replace(".ali", "")
        mod_name = str(basename) + '-mod.ali'
        mod_file = open(mod_name, "w")
        alifile = open(ali)
        for line in alifile:
            record = line.replace("\n", "")
            if '>' in line:
                defline = record.replace(">", "")
                split_list = defline.split("@")
                long = split_list[0]
                ID = split_list[1]
                short = idm_ofs[long]
                mod_file.write('>' + str(short) + '@' + str(ID) + "\n")
            else:
                mod_file.write(str(record) + "\n")


if __name__ == '__main__':
    main()
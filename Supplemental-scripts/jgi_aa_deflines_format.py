#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""jgi_aa_deflines_format.py"""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2021, Sciensano"
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

    jgi_file = open(main_file)
    basename = main_file.replace(".faa", "")
    mod_name = str(basename) + '-format.faa'
    mod_file = open(mod_name, "w")
    for line in jgi_file:
        record = line.replace("\n", "")
        if '>' in line:
            defline = record.replace(">", "")
            split_list = defline.split("|")
            genome = split_list[1]
            ID = split_list[2]
            string = str(genome) + '-' + str(ID)
            mod_file.write('>' + str(string) + "\n")
        else:
            mod_file.write(str(record) + "\n")


if __name__ == '__main__':
    main()
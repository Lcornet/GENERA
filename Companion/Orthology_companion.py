#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Orthology companiion."""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2021, University of Liège"
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
#taxid
@click.option('--duplication', default='0')
@click.option('--unwanted', default='0')
@click.option('--duplicationpercent', default='200')
@click.option('--presence', default='100')

def main(main_file, duplication, presence, duplicationpercent, unwanted):

    #open list of organism
    orgs_of = {}

    infile = open(main_file)
    for line in infile:
        record = line.replace("\n", "")
        orgs_of[record] = 1
    orgs_number = len(orgs_of)

    #open OGs to define core gene
    OG_list = glob.glob("OG*.fa")
    out_file = open('core-OG.list', "w")
    for OG in OG_list:
        #define
        count = 0
        unwanted_number = 0
        core_of = {}
        seen_of = {}

        OGfile = open(OG)
        for line in OGfile:
            file_record = line.replace("\n", "")
            if '>' in line:
                defline = file_record.replace(">", "")
                split_list = defline.split("|")
                defline = split_list[0]
                #load
                if (defline in seen_of):
                    seen_of[defline] += 1
                else:
                    seen_of[defline] = 1
                count += 1
                #is it in orgs for core genes ?
                if (defline in orgs_of):
                    if (defline in core_of):
                        core_of[defline] += 1
                    else:
                        core_of[defline] = 1
                else:
                    unwanted_number += 1

        #compute stat for this OG
        dups = 0
        for defline in seen_of:
            number = seen_of[defline]
            if (number > 1):
                dups += 1
        dups_percentage = (float(dups) / float(count)) * 100
        core = len(core_of)
        orgs_percentage = (float(core) / float(orgs_number)) * 100

        #print OG
        dupli = int(duplication)
        dupliPercent = float(duplicationpercent)
        pres = float(presence)
        
        if (orgs_percentage >= pres):
            #pass presence check
            if (unwanted_number <= int(unwanted)):
                #less unwated or than desired
                if (dups <= dupli):
                    #pass duplication absolute filter
                    if (dups_percentage <= dupliPercent):
                        out_file.write( str(OG) + "\n")

if __name__ == '__main__':
    main()
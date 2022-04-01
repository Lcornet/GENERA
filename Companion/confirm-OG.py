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
import glob

@click.command()
###ARGUMENT#####
#main file - depending of the mode can be checkm result or rnammer output
@click.argument('main_file', type=click.Path(exists=True,readable=True))

def main(main_file):


    #Openlist file
    list_of = {}
    list = open(main_file)
    for line in list:
        record = line.replace("\n", "")
        list_of[record]=1

    #open final list
    list_out = open('verified-specific.list', "w")

    #Open OG
    OG_list = glob.glob("*OG*-GENERA.faa")
    for OG in OG_list:
        seen_genome = {} 
        unwanted_genome = {} 
        printtag = 0
        OGfile = open(OG)
        #open sub fro writing
        subname = OG.replace("GENERA.faa", "specific.faa")
        file_out = open(subname, "w")
        for line in OGfile:
            record = line.replace("\n", "")
            if '>' in line:
                printtag = 0
                newtag = 0
                if '#NEW#' in line:
                    newtag = 1
                defline = record.replace(">", "")
                split_list = defline.split("@")
                genome = split_list[0]
                genome = genome.replace(" ", "_")
                id = split_list[1]
                #Print the sub sampled gene
                if genome in list_of:
                    #Genome required for specific content
                    if (newtag == 0):
                        #Genome alredy present; printed
                        file_out.write('>' + str(genome) + '@' + str(id) + "\n")
                        seen_genome[genome] = 1
                        printtag = 1
                    else:
                        #Genome added by forty-tow, check that it is not already printed
                        if (genome not in seen_genome):
                            file_out.write('>' + str(genome) + '@' + str(id) + "\n")
                            seen_genome[genome] = 1
                            printtag = 1
                else:
                    unwanted_genome[genome] = 1
            else:
                if (printtag == 1):
                    file_out.write(str(record) + "\n")
        #Add to the final specific list ?
        requiredlen = len(list_of)
        fastasize = len(seen_genome)
        unwantedsize = len(unwanted_genome)
        if ((fastasize == requiredlen) and (unwantedsize == 0)):
            list_out.write(str(subname) + "\n")

if __name__ == '__main__':
    main()
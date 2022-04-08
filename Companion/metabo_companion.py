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
from collections import defaultdict

@click.command()
###ARGUMENT#####
#main file - depending of the mode can be checkm result or rnammer output
@click.argument('main_file', type=click.Path(exists=True,readable=True))
#option
@click.option('--mode', default='no', help='mode')

def main(main_file, mode):

    #Open anvio files
    anvio_list = glob.glob("*metabo_modules.txt")
    #declare
    metabo_of = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
    seen_module_of = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
    seen_genome_of = {}

    for anvio in anvio_list:
        anviofile = open(anvio)
        for line in anviofile:
            record = line.replace("\n", "")
            if ('unique_id' in record):
                continue
            else:
                split_list = record.split("\t")
                ID = split_list[0]
                genome = split_list[1]
                kegg_module = split_list[2]
                module_name = split_list[3]
                chunks = module_name.split(",")
                short_module_name = chunks[0]
                short_module_name = short_module_name.replace(" ", "-")
                module_category = split_list[5]
                module_category = module_category.replace(" ", "_")
                module_subcategory = split_list[6]
                module_definition = split_list[7]
                module_completeness = split_list[8]
                module_complete = split_list[9]

                #load
                kegg_ID = str(kegg_module) + '_' + str(short_module_name)
                #First, no-zero mode
                if (mode == 'no-zero'):
                    if (module_complete == 'True'):
                        metabo_of[module_category][genome][kegg_ID] = 1
                        seen_genome_of[genome]= 1
                        seen_module_of[module_category][kegg_ID] = 1
                elif (mode == 'show-zero'):
                    seen_genome_of[genome]=1
                    seen_module_of[module_category][kegg_ID] = 1
                    if (module_complete == 'True'):
                        metabo_of[module_category][genome][kegg_ID] = 1
                    elif (module_complete == 'False'):
                        metabo_of[module_category][genome][kegg_ID] = 0

    #print
    for module_category in metabo_of:
        #open outfile for category
        out = str(module_category) + '-infile.txt'
        out_file = open(out, "w")
        out_file.write('KEGG' + "\t" + "GENOME" + "\t" + "PRESENCE" + "\n")
        #denest
        category_of = metabo_of[module_category]
        module_of = seen_module_of[module_category]
        #check the presence of the module in the genome
        for genome in seen_genome_of:
            genome_of = category_of[genome]
            #Lopp in seen kegg:
            for kegg_ID in module_of:
                if (kegg_ID in genome_of):
                    num = genome_of[kegg_ID]
                    if (num == 1):
                        out_file.write(str(kegg_ID) + "\t" + str(genome) + "\t" + '1' + "\n")
                    elif (num == 0):
                        out_file.write(str(kegg_ID) + "\t" + str(genome) + "\t" + '0' + "\n")
                else:
                    out_file.write(str(kegg_ID) + "\t" + str(genome) + "\t" + '0' + "\n")


if __name__ == '__main__':
    main()
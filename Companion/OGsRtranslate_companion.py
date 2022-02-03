#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Companion"""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2020, University of Liège"
__version__ = "1.0.0"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"

from fileinput import filename
import click

@click.command()
###ARGUMENT#####
#main file - depending of the mode can be checkm result or rnammer output
@click.argument('main_file', type=click.Path(exists=True,readable=True))
#taxid
@click.option('--mode', default='no', help='mode')

def main(main_file, mode):
		
	#read fasta file
    if (mode == 'yaml'):
        file_out = open("yaml.part2", "w")
        #Open file matrix
        infile = open(main_file)
        for line in infile:
            sum_record = line.replace("\n", "")
            split_list = sum_record.split(" ")
            filename = split_list[0]
            file_out.write('  - org: Genus species_' + str(filename) +  "\n")
            file_out.write('    banks:' + "\n")
            #file_out.write('        - ' +  str(filename) + '.fna' + "\n")  
            file_out.write('        - ' +  str(filename) + "\n") 
    
    if (mode == 'fsp'):
        chunks = main_file.split(".")
        basename = chunks[0]
        name = str(basename) + '-sp.ali'
        file_out = open(name, "w")
        infile = open(main_file)
        for line in infile:
            sum_record = line.replace("\n", "")
            if ('>' in sum_record):
                #Get GCA
                defline = sum_record.replace(">", "")
                split_list = defline.split("@")
                GCA = split_list[0]
                GCA = GCA.replace(" ", "_")
                ID = split_list[1]
                file_out.write('>Genus species_' + str(GCA) + '@' + str(ID) + "\n")
            else:
                file_out.write(str(sum_record) + "\n")
    
    if (mode == 'restore'):
        chunks = main_file.split(".")
        basename = chunks[0]
        basename = basename.replace("-sp", "")
        name = str(basename) + '.ali'
        file_out = open(name, "w")
        infile = open(main_file)
        for line in infile:
            sum_record = line.replace("\n", "")
            if ('>' in sum_record):
                #Get GCA
                defline = sum_record.replace(">Genus species_", "")
                split_list = defline.split("@")
                GCA = split_list[0]
                ID = split_list[1]
                file_out.write('>' + str(GCA) + '@' + str(ID) + "\n")
            else:
                file_out.write(str(sum_record) + "\n")




if __name__ == '__main__':
    main()
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Kraken-ids-creation.py: Creation of ids to add genome in kraken db library."""
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
#taxid
@click.option('--mode', default='no', help='mode')

def main(main_file, mode):
		
	#read fasta file
    if (mode == 'third'):
        loop = 0

        #Open file matrix
        infile = open(main_file)
        for line in infile:
            sum_record = line.replace("\n", "")
            loop += 1
            split_list = sum_record.split(" ")
            length = split_list[1]
            if (loop == 1):
                #Print partition file
                third_out = open("partition.txt", "w")
                third_out.write('DNA, p1=1-' + str(length) + '\\3,2-' + str(length) + '\\3' + "\n")
                third_out.write('DNA, p2=3-' + str(length) + '\\3' + "\n")
                break
    
    if (mode == 'two'):

        length = 'NA'
        loop = 0

        #Open file matrix
        infile = open(main_file)
        for line in infile:
            sum_record = line.replace("\n", "")
            loop += 1
            split_list = sum_record.split(" ")
            length = split_list[1]
            if (loop == 1):
                break

        #Create blocks file
        two_out = open("blocks", "w")
        progress = 0
        loop = 0
        while progress < int(length):
            progress += 1
            loop += 1
            if (loop == 3):
                loop = 0
                two_out.write("\n")
            else:
                two_out.write(str(progress) + " ")

if __name__ == '__main__':
    main()
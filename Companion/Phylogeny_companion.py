#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Kraken-ids-creation.py: Creation of ids to add genome in kraken db library."""
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
    
    if (mode == 'IDM'):

        #IDL part
        infile = open(main_file)
        idl_of = {}
        for line in infile:
            record = line.replace("\n", "")
            if ('#' in record):
                continue
            else:
                idl_of[record] = 1

        #IDM part
        infile = open('IDM')
        temp_idm = open("idm.temp", "w")
        temp_idm.write("#" + "\n")
        temp_idm.write("#" + "\n")

        for line in infile:
            record = line.replace("\n", "")
            if ('#' in record):
                continue
            else:
                split_list = record.split("\t")
                long = split_list[0]
                short = split_list[1]
                newshort = 'NA'
                for id in idl_of:
                    chunks = id.split("@")
                    shortID = chunks[0]
                    if (shortID == short):
                        newshort = id
                temp_idm.write(str(long) + "\t" + str(newshort) + "\n")
    
    if (mode == 'iter'):

        ali_list = glob.glob("*.ali")
        for ali in ali_list:
            alifile = open(ali)
            aliiter = ali.replace(".ali", "-iter.ali")
            ali_iter = open(aliiter, "w")
            num = 1
            for line in alifile:
                record = line.replace("\n", "")
                if ('>') in record:
                    seq = str(record) + '@' + str(num)
                    num += 1
                    ali_iter.write(str(seq) + "\n")
                else:
                    ali_iter.write(str(record) + "\n")

if __name__ == '__main__':
    main()
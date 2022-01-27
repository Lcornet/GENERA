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
@click.option('--taxa', default='no', help='Choice between four taxa levels: phylum, class, order, family')
@click.option('--refgroup', default='no', help='Name of the group')

def main(main_file, mode, taxa, refgroup):
		
	#read fasta file
    if (mode == 'sum'):
        #define if it's from genbank or refseq
        origin = main_file
        origin = origin.replace("_sum.txt", "")

        if (origin == "refseq"):
            #outfile
            sum_out = open("refseq_sum-filt.txt", "w")

            #Open checkm result
            infile = open(main_file)

            for line in infile:
                sum_record = line.replace("\n", "")
                if ('#' in sum_record):
                    continue
                else:
                    split_list = sum_record.split("\t")
                    link = split_list[19]
                    if (link != 'na'):
                        sum_out.write(sum_record + "\n")
        
        elif (origin == 'genbank'):
            #outfile
            sum_out = open("genbank_sum-filt.txt", "w")

            #Open checkm result
            infile = open(main_file)

            for line in infile:
                sum_record = line.replace("\n", "")
                if ('#' in sum_record):
                    continue
                else:
                    if ('http' in sum_record): 
                        #print(sum_record)
                        split_list = sum_record.split("\t")
                        link = split_list[19]
                        if (link != 'na'):
                            sum_out.write(sum_record + "\n")
                    else:
                        continue
    
    elif (mode == 'fetch'):
        #define if it's from genbank or refseq
        origin = main_file
        origin = origin.replace(".tax", "")

        if (origin == "GCF"):
            #Open fetch tax file
            infile = open(main_file)

            GCFlist_of = {}
            for line in infile:
                fetch_record = line.replace("\n", "")
                split_list = fetch_record.split("\t")
                GC = split_list[0]
                levels = split_list[3]
                #keep only GC ids without number
                split_list = GC.split(".")
                GC = split_list[0]
                #separate taxa taxonomic levels
                if (';' not in levels):
                    continue
                else: 
                    split_list = levels.split(";")
                    phylum = split_list[0]
                    phylum = phylum.replace(" ", "")
                    clas = split_list[1]
                    clas = clas.replace(" ", "")
                    order = split_list[2]
                    order = order.replace(" ", "")
                    family = split_list[3]
                    family = family.replace(" ", "")
                    genus = split_list[4]
                    genus = genus.replace(" ", "")

                    if (taxa == 'phylum'):
                        group = phylum
                        if (refgroup in group):
                            GCFlist_of[GC] = 1
                    elif (taxa == 'class'):
                        group = clas
                        if (refgroup in group):
                            GCFlist_of[GC] = 1
                    elif (taxa == 'order'):
                        group = order
                        if (refgroup in group):
                            GCFlist_of[GC] = 1
                    elif (taxa == 'family'):
                        group = family
                        if (refgroup in group):
                            GCFlist_of[GC] = 1
                    elif (taxa == 'genus'):
                        group = genus
                        if (refgroup in group):
                            GCFlist_of[GC] = 1
                    
            
            #print GC list
            fetch_out = open("GCF.refgroup.uniq", "w")
            for GC in GCFlist_of:
                fetch_out.write(str(GC) + "\n")
        
        elif (origin == "GCA"):
            #open refgroupRefseqGC to lod ids of refseq genomes already in the db
            refseqList = open('GCF.refgroup.uniq')

            GCFlist_of = {}
            for line in refseqList:
                reflist_record = line.replace("\n", "")
                ID = reflist_record.replace("GCF_", "")
                GCFlist_of[ID] = 1

            #Open fetch tax file
            infile = open(main_file)

            GCAlist_of = {}
            for line in infile:
                fetch_record = line.replace("\n", "")
                split_list = fetch_record.split("\t")
                GC = split_list[0]
                levels = split_list[3]
                #keep only GC ids without number
                split_list = GC.split(".")
                GC = split_list[0]
                ID = GC.replace("GCA_", "")
                #separate taxa taxonomic levels
                if (';' not in levels):
                    continue
                else: 
                    split_list = levels.split(";")
                    phylum = split_list[0]
                    phylum = phylum.replace(" ", "")
                    clas = split_list[1]
                    clas = clas.replace(" ", "")
                    order = split_list[2]
                    order = order.replace(" ", "")
                    family = split_list[3]
                    family = family.replace(" ", "")
                    genus = split_list[4]
                    genus = genus.replace(" ", "")

                    if (taxa == 'phylum'):
                        group = phylum
                        if ((refgroup in group) and (ID not in GCFlist_of)):
                            GCAlist_of[GC] = 1
                    elif (taxa == 'class'):
                        group = clas
                        if ((refgroup in group) and (ID not in GCFlist_of)):
                            GCAlist_of[GC] = 1
                    elif (taxa == 'order'):
                        group = order
                        if ((refgroup in group) and (ID not in GCFlist_of)):
                            GCAlist_of[GC] = 1
                    elif (taxa == 'family'):
                        group = family
                        if ((refgroup in group) and (ID not in GCFlist_of)):
                            GCAlist_of[GC] = 1
                    elif (taxa == 'genus'):
                        group = genus
                        if ((refgroup in group) and (ID not in GCFlist_of)):
                            GCAlist_of[GC] = 1
    
            #print GC list
            fetch_out = open("GCA.refgroup.uniq", "w")
            for GC in GCAlist_of:
                fetch_out.write(str(GC) + "\n")

if __name__ == '__main__':
    main()
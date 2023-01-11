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
@click.option('--mode', default='manytomany', help='mode')
@click.option('--submode', default='dist', help='submode')
@click.option('--anitool', default='fastANI', help='anitool')

def main(main_file, mode, submode, anitool):

    if ((mode == 'manytomany') or (mode == 'onetomany')):

        #Open main file = ANI values
        ANI_of = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        PROX_of = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        DUP_of = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        anifile = open(main_file)
        for line in anifile:
            record = line.replace("\n", "")
            split_list = record.split("\t")
            genome1 = split_list[0]
            genome1 = genome1.replace(".fna", "")
            genome2 = split_list[1]
            genome2 = genome2.replace(".fna", "")
            ANI = split_list[2]
            ANI_of[genome1][genome2] = ANI
            #When identical genomes, tell user
            if (genome1 in ANI_of):
                denest_of = PROX_of[genome1]
                if (ANI in denest_of):
                    #Duplicate hits in ANI
                    duplicate = denest_of[ANI]
                    DUP_of[genome1][genome2] = duplicate
                PROX_of[genome1][ANI] = genome2
            else:
                PROX_of[genome1][ANI] = genome2

        #Open list of genome
        genome_lists=[]
        listfile = open('list')
        for line in listfile:
            record = line.replace("\n", "")
            genome_lists.append(record)

        #Print matrix
        out_file = open('ANI-matrix.txt', "w")
        #Master loop
        for master in genome_lists:
            out_file.write(str(master) + "\t")
            denest_of = ANI_of[master]
            #Slave loop
            for slave in genome_lists:
                if (slave in denest_of):
                    ANI = denest_of[slave]
                else:
                    ANI = 'NA'
                out_file.write(str(ANI) + "\t")
            #Close the line
            out_file.write("\n")

        #Print dist matrix
        dist_file = open('ANI-dist-matrix.txt', "w")
        #Master loop
        distlen = len(genome_lists)
        dist_file.write(str(distlen) + "\n")
        for master in genome_lists:
            dist_file.write(str(master) + "\t")
            denest_of = ANI_of[master]
            #Slave loop
            listlen = len(genome_lists)
            listcount = 0
            for slave in genome_lists:
                listcount += 1
                if (slave in denest_of):
                    ANI = denest_of[slave]
                    if (anitool == 'fastANI'):
                        ANI = 1 - (float(ANI) / 100)
                    elif (anitool == 'orthoani'):
                        ANI = 1 - (float(ANI))
                else:
                    ANI = 'NA'
                #Don't print tabluation for last element
                if (listcount == listlen):
                    dist_file.write(str(ANI))
                else:
                    dist_file.write(str(ANI) + "\t")
            #Close the line
            dist_file.write("\n")
        
        #Print ggplot2 infile
        gg_file = open('ANI-infile.txt', "w")
        gg_file.write('GENOME1' + "\t" + 'GENOME2' + "\t" + 'ANI' + "\n")
        #Master loop
        for master in genome_lists:
            denest_of = ANI_of[master]
            #Slave loop
            for slave in genome_lists:
                gg_file.write(str(master) + "\t")
                gg_file.write(str(slave) + "\t")
                if (slave in denest_of):
                    ANI = denest_of[slave]
                else:
                    ANI = 'NA'
                gg_file.write(str(ANI) + "\n")
    
    #Classical moden all or shortlist mode
    if (mode == 'manytomany'):
        #Print proximity
        #Master loop
        for master in genome_lists:
            cat = str(master) + '-closest-ANI.list'
            ANI_file = open(cat, "w")
            #First print any identical genomes
            if (master in DUP_of):
                denest_of = DUP_of[master]
                for duplicate in denest_of:
                    identical = denest_of[duplicate]
                    ANI_file.write(str(identical) + ' seems to be identical to ' + str(duplicate) + "\n" )

            #closet part
            denest_of = PROX_of[master]
            #Sort denest of master based on ANI
            for ani in denest_of:
                print(ani)
                slave = denest_of[ani]
                ANI_file.write(str(slave) + "\t" + str(ani) + "\n")

    elif (mode == 'onetomany'):

        #Open list of genome
        shortgenome_lists=[]
        listfile = open('shortlist')
        for line in listfile:
            record = line.replace("\n", "")
            shortgenome_lists.append(record)

        #Master loop
        for master in shortgenome_lists:
            cat = str(master) + '-closest-ANI.list'
            ANI_file = open(cat, "w")
            #First print any identical genomes
            if (master in DUP_of):
                denest_of = DUP_of[master]
                for duplicate in denest_of:
                    identical = denest_of[duplicate]
                    ANI_file.write(str(identical) + ' seems to be identical to ' + str(duplicate) + "\n" )

            #closet part
            denest_of = PROX_of[master]
            #Sort denest of master based on ANI
            for ani in denest_of:
                #print(ani)
                slave = denest_of[ani]
                ANI_file.write(str(slave) + "\t" + str(ani) + "\n")
    
    elif (mode == 'idm'):

        #Read idm
        idmfile = open('file.idm')
        idm_of = {}
        for line in idmfile:
            record = line.replace("\n", "")
            if ('#' in record):
                continue
            else:
                split_list = record.split("\t")
                long = split_list[0]
                short = split_list[1]
                idm_of[short] = long
            
        if (submode == 'dist'):
            #Change the name in newick file based on IDM
            distfile = open('ANI-dist-matrix.txt')
            PHY_file = open('ANI.phylip', "w")
            for line in distfile:
                record = line.replace("\n", "")
                split_list = record.split("\t")
                lenlist = len(split_list)
                countlist = 0
                for f in split_list:
                    field = 'NA'
                    if f in idm_of:
                        field = idm_of[f]
                    else:
                        field = f
                    #write part
                    countlist += 1
                    if (countlist == lenlist):
                        PHY_file.write(str(field) + "\n")
                    else:
                        PHY_file.write(str(field) + "\t")

        elif (submode == 'heatmap'):
            #Change the name in newick file based on IDM
            heatfile = open('ANI-infile.txt')
            HEAT_file = open('ANI-heat.txt', "w")
            for line in heatfile:
                record = line.replace("\n", "")
                split_list = record.split("\t")
                lenlist = len(split_list)
                countlist = 0
                for f in split_list:
                    field = 'NA'
                    if f in idm_of:
                        field = idm_of[f]
                    else:
                        field = f
                    #write part
                    countlist += 1
                    if (countlist == lenlist):
                        HEAT_file.write(str(field) + "\n")
                    else:
                        HEAT_file.write(str(field) + "\t")            


if __name__ == '__main__':
    main()
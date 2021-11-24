#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Assembly companiion."""
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
@click.option('--mode', default='no', help='mode')

def main(main_file, mode):

    #Compute binning percentage

    #open main file
    length_of = {}
    defline_number = 0
    tot_len = 0
    count = 0
    def_length = 0
    defline = 'NA'

    infile = open(main_file)
    for line in infile:
        fasta_record = line.replace("\n", "")
        if ('>' in fasta_record):
            if (count >= 1):
                #load previous defline
                length_of[defline] = def_length
                #reset
                def_length = 0
                defline = 'NA'
            defline = fasta_record.replace(">", "")
            defline_number += 1
            count += 1
        else:
            length = len(fasta_record)
            def_length = def_length + length
            tot_len += length
    #last contig
    length_of[defline] = def_length

    
    #Open bins
    if (mode == 'all'):
        #metabat
        metabat_of = {}
        metabat_list = glob.glob("METABAT_bin*.fa")
        for metabat in metabat_list:
            Mfile = open(metabat)
            for line in Mfile:
                file_record = line.replace("\n", "")
                if '>' in line:
                    seq = file_record.replace(">", "")
                    metabat_of[seq] = 1

        #concoct
        concoct_of = {}
        concoct_list = glob.glob("CONCOCT_bin*.fa")
        for concoct in concoct_list:
            Cfile = open(concoct)
            for line in Cfile:
                file_record = line.replace("\n", "")
                if '>' in line:
                    seq = file_record.replace(">", "")
                    concoct_of[seq] = 1

        #compute percentage
        binned_contigs = 0
        binned_length = 0
        for contig in length_of:
            if (contig in metabat_of):
                binned_contigs += 1
                length = length_of[contig]
                binned_length += length
        #stats
        metabat_totContigs = (float(binned_contigs) / float(defline_number)) * 100
        metabat_totlen = (float(binned_length) / float(tot_len)) * 100

        #compute percentage
        binned_contigs = 0
        binned_length = 0
        for contig in length_of:
            if (contig in concoct_of):
                binned_contigs += 1
                length = length_of[contig]
                binned_length += length
        #stats
        concoct_totContigs = (float(binned_contigs) / float(defline_number)) * 100
        concoct_totlen =  (float(binned_length) / float(tot_len)) * 100

        out_file = open('binned.info', "w")
        #open GENERA-Assembler.log
        info_file = open('GENERA-Assembler.log')
        for line in info_file:
            file_record = line.replace("\n", "")
            out_file.write(str(file_record) + "\n")
        out_file.write('GENERA info: Metabat binned ' + str(metabat_totContigs) + ' % of contigs and ' + str(metabat_totlen) + ' of the length of the metagenome' + "\n")
        out_file.write('GENERA info: CONCOCT binned ' + str(concoct_totContigs) + ' % of contigs and ' + str(concoct_totlen) + ' of the length of the metagenome' + "\n")
    
    elif (mode == 'metabat'):
        #metabat
        metabat_of = {}
        metabat_list = glob.glob("METABAT_bin*.fa")
        for metabat in metabat_list:
            Mfile = open(metabat)
            for line in Mfile:
                file_record = line.replace("\n", "")
                if '>' in line:
                    seq = file_record.replace(">", "")
                    metabat_of[seq] = 1
        
        #compute percentage
        binned_contigs = 0
        binned_length = 0
        for contig in length_of:
            if (contig in metabat_of):
                binned_contigs += 1
                length = length_of[contig]
                binned_length += length

        #stats
        metabat_totContigs = (float(binned_contigs) / float(defline_number)) * 100
        metabat_totlen = (float(binned_length) / float(tot_len)) * 100

        out_file = open('binned.info', "w")
        #open GENERA-Assembler.log
        info_file = open('GENERA-Assembler.log')
        for line in info_file:
            file_record = line.replace("\n", "")
            out_file.write(str(file_record) + "\n")
        out_file.write('GENERA info: Metabat binned ' + str(metabat_totContigs) + ' % of contigs and ' + str(metabat_totlen) + ' % of the length of the metagenome' + "\n")


    elif (mode == 'concoct'):
        #concoct
        concoct_of = {}
        concoct_list = glob.glob("CONCOCT_bin*.fa")
        for concoct in concoct_list:
            Cfile = open(concoct)
            for line in Cfile:
                file_record = line.replace("\n", "")
                if '>' in line:
                    seq = file_record.replace(">", "")
                    concoct_of[seq] = 1
        
        #compute percentage
        binned_contigs = 0
        binned_length = 0
        for contig in length_of:
            if (contig in concoct_of):
                binned_contigs += 1
                length = length_of[contig]
                binned_length += length
        #stats
        concoct_totContigs = (float(binned_contigs) / float(defline_number)) * 100
        concoct_totlen = (float(binned_length) / float(tot_len)) * 100       

        out_file = open('binned.info', "w")
        #open GENERA-Assembler.log
        info_file = open('GENERA-Assembler.log')
        for line in info_file:
            file_record = line.replace("\n", "")
            out_file.write(str(file_record) + "\n")
        out_file.write('GENERA info: CONCOCT binned ' + str(concoct_totContigs) + ' % of contigs and ' + str(concoct_totlen) + ' % of the length of the metagenome' + "\n")
    
    

if __name__ == '__main__':
    main()
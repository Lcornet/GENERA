#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""change-ext.py: Change ext of fasta/ali files."""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2020, University of Liège"
__version__ = "1.0.1"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"


def main():

    #tool = 'Checkm_contamination'

    levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    tools_of = ['CheckM', 'BUSCO', 'GUNC', 'Physeter', 'Kraken2', 'CheckM2']


    #Open out wiles
    conta_file = open('conta.txt', "w")
    conta_file.write('level' + "\t" + 'percent' + "\t" + 'ID' + "\n")
    #results_file = open('results.txt', "w")
    #results_file.write('level' + "\t" + 'percent' + "\n")

    #open chime
    for level in levels:
        name = 'chimeric-genomes_' + str(level) + '.list'
        levelfile = open(name)
        for line in levelfile:
            if ('available' in line):
                continue
            record = line.replace("\n", "")
            split_list = record.split("\t")
            percent = split_list[12]
            percent = percent.replace("Chimeric_level=", "")
            ID = split_list[4]
            ID = ID.replace("creation of chimeric genome: ID=", "")
            conta_file.write(str(level) + "\t" + str(percent) + "\t" + str(ID) + "\n")

    #open results
    for tool in tools_of:
        #Open outfile
        cat = 'results' + str(tool) + '.txt'
        results_file = open(cat, "w")
        results_file.write('level' + "\t" + 'percent' + "\t" + 'category' + "\t" + 'ID' + "\n")

        for level in levels:
            #name = 'chimeric-genomes_' + str(level) + '.list'
            name = 'GENERA-contams-' + str(level) + '.table'
            resultfile = open(name)
            for line in resultfile:
                record = line.replace("\n", "")
                if ('Genome' in record):
                    continue
                split_list = record.split("\t")
                ID = split_list[0]
                ID = ID.replace("chimeric-","")
                ID = ID.replace("-abbr","")
                Checkm_contamination = split_list[2]
                Checkm_str = split_list[3]
                Busco_duplicate = split_list[6]
                GUNC_CSS = float(split_list[7]) * 100
                GUNC_RRS = float(split_list[8]) * 100
                GUNC_conta = float(split_list[10]) * 100
                Physeter_contamination = split_list[12]        
                Kraken_contamination = split_list[14]
                Checkm2_contamination = split_list[16]

                #Pass in each tool and print
                if (tool == 'CheckM'):
                    results_file.write(str(level) + "\t" + str(Checkm_contamination) + "\t" + 'CK-conta' + "\t" + str(ID) + "\n")
                    results_file.write(str(level) + "\t" + str(Checkm_str) + "\t" + 'CK-strh' + "\t" + str(ID) + "\n")
                elif (tool == 'BUSCO'):
                    results_file.write(str(level) + "\t" + str(Busco_duplicate) + "\t" + 'BUSCO-dup' + "\t" + str(ID) + "\n")
                elif (tool == 'GUNC'):
                    #results_file.write(str(level) + "\t" + str(GUNC_CSS) + "\t" + 'GUNC-CSS' + "\t" + str(ID) + "\n")
                    #results_file.write(str(level) + "\t" + str(GUNC_RRS) + "\t" + 'GUNC-RRS' + "\t" + str(ID) + "\n")
                    results_file.write(str(level) + "\t" + str(GUNC_conta) + "\t" + 'GUNC-conta' + "\t" + str(ID) + "\n")
                elif (tool == 'Physeter'):
                    results_file.write(str(level) + "\t" + str(Physeter_contamination) + "\t" + 'Physeter-conta' + "\t" + str(ID) + "\n")
                elif (tool == 'Kraken2'):
                    results_file.write(str(level) + "\t" + str(Kraken_contamination) + "\t" + 'Kraken2-conta' + "\t" + str(ID) + "\n")
                elif (tool == 'CheckM2'):
                    results_file.write(str(level) + "\t" + str(Checkm2_contamination) + "\t" + 'CheckM2-conta' + "\t" + str(ID) + "\n")


if __name__ == '__main__':
    main()
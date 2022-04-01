#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""anvio-to-BMC: modify deflines to match BMC"""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2021, Sciensano"
__version__ = "1.0.0"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"


import click
import glob

@click.command()
###ARGUMENT#####
#main file - Postitive List of Orgs 
#taxid
@click.option('--idm', default='no', help='path-to-idm')

def main(idm):

    #Read the list or orf and store
    idmfile = open(idm)
    idm_of = {}
    for line in idmfile:
        list_record = line.replace("\n", "")
        split_list = list_record.split("\t")
        long = split_list[0]
        abbr = split_list[1]
        idm_of[abbr] = long

    #print(idm_of)


    OG_list = glob.glob("*OG.fasta")
    for OG in OG_list:
        Ortho = OG
        Ortho = Ortho.replace("_OG.fasta", "")
        bmcname = str(Ortho) + '-BMC_OG.fasta'
        out_file = open(bmcname, "w")
        OGfile = open(OG)
        for line in OGfile:
            file_record = line.replace("\n", "")
            if '>' in line:
                split_list = file_record.split("|")
                chunk = split_list[2]
                id = split_list[0]
                id = id.replace(">", "")
                split_list = chunk.split(":")
                org = split_list[1]
                string = 'na'
                if (idm == 'no'):
                    string = '>' + str(org) + '@' + str(id)
                else:
                    org = idm_of[org] #switch to full name
                    string = '>' + str(org) + '@' + str(id)
                out_file.write(str(string) + "\n")
            else:
                out_file.write(str(file_record) + "\n")

if __name__ == '__main__':
    main()

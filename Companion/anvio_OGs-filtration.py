#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""anvio_OGs-filtration: Create a list of OGs"""
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
@click.argument('main_file', type=click.Path(exists=True,readable=True))
#taxid
@click.option('--mode', default='anvio', help='OG inference mode: anvio or OF')
@click.option('--pfilter', default='no', help='presence filter: yes or no')
@click.option('--fraction', default='1', help='fraction of of orgs from the positive list needed to conserve a file')
@click.option('--unwanted', default='0', help='Maximum of unwanted org in a file to be conserved')
@click.option('--cfilter', default='no', help='copy filter: yes or no')
@click.option('--maxcopy', default='1', help='Mean of maximum of sequence of the same org in a file to be conserved')
@click.option('--hfilter', default='no', help='homogeneity filter (anvio only): yes or no')
@click.option('--hindex', default='homogeneity.txt', help='homogeneity index (anvio only): path')
@click.option('--maxfunctionalindex', default='0.95', help='Maximum functional index of anvio (1 perfect) for a file to be conserved')
@click.option('--maxgeometricindex', default='0.95', help='Maximum geometric index of anvio (1 perfect) for a file to be conserved')
@click.option('--cogtable', default='../SUMMARY/LMG_Pan_gene_clusters_summary.txt', help='path to cog table')

def main(main_file, mode, pfilter, fraction, unwanted, cfilter, maxcopy, hfilter, hindex, maxfunctionalindex, maxgeometricindex, cogtable):

    #open outfile
    out_file = open('filtered-OG.list', "w")

    #Read the list or orf and store
    orgs_of = {}
    listfile = open(main_file)
    for line in listfile:
        list_record = line.replace("\n", "")
        split_list = list_record.split("\t")
        org = split_list[0]
        orgs_of[org] = 1
    orgNumber = len(orgs_of)
    #print(orgs_of)
    
    #Get list of OGs
    OG_list = glob.glob("*OG.fasta")

    #Filter OGs based on positive list and fraction
    if (pfilter == 'yes'):
        tempOGlist_of = []
        for OG in OG_list:
            orgCount = 0
            counted_orgs = []
            seen_orgs = []
            unwantedorgCount = 0
            OGfile = open(OG)
            for line in OGfile:
                file_record = line.replace("\n", "")
                #deflines
                if '>' in line:
                    org = 'na'
                    if (mode == 'anvio'):
                        split_list = file_record.split("|")
                        chunk = split_list[2]
                        split_list = chunk.split(":")
                        org = split_list[1]
                    elif (mode == 'OG'):
                        split_list = file_record.split("|")
                        org = split_list[0]
                        org = org.replace(">", "")
                    #print(org)
                    if ((org in orgs_of) and (org not in counted_orgs)): #Counte only one time an organism
                        orgCount += 1
                        counted_orgs.append(org)
                    else:
                        if (org not in seen_orgs):
                            unwantedorgCount += 1
                            seen_orgs.append(org)
            #End of OG, compute presence fraction
            presenceFraction = 0
            if (orgCount >= 1): #illegal division by zero
                presenceFraction = float(orgCount)/float(orgNumber)
            if ((presenceFraction >= float(fraction)) and (unwantedorgCount <= int(unwanted))):
                tempOGlist_of.append(OG)
        #Replace OG list with the new one
        totOGcount = len(OG_list)
        OG_list = tempOGlist_of
        keepOGcount = len(OG_list)
        lostOG = int(totOGcount) - int(keepOGcount)
        #Print stat
        string1= 'pFilter activated: Conserve only OG with ' + str(fraction) + ' X the tot number of positive org and ' + str(unwanted) + ' unwated org.'
        out_file.write('#' + str(string1) + "\n")
        string2= str(keepOGcount) + ' conserved OGs and ' + str(lostOG) + ' lost OGs.'
        out_file.write('#' + str(string2) + "\n")
        print(string1)
        print(string2)
    
    #Filter based on copy number
    if (cfilter == 'yes'):  
        tempOGlist_of = []
        for OG in OG_list:
            seen_orgs = {}
            OGfile = open(OG)
            for line in OGfile:
                file_record = line.replace("\n", "")
                #deflines
                if '>' in line:
                    org = 'na'
                    if (mode == 'anvio'):
                        split_list = file_record.split("|")
                        chunk = split_list[2]
                        split_list = chunk.split(":")
                        org = split_list[1]
                    elif (mode == 'OG'):
                        split_list = file_record.split("|")
                        org = split_list[0]
                        org = org.replace(">", "")
                    #load org in dict
                    if (org not in seen_orgs):
                        seen_orgs[org] = 1
                    else:
                        seen_orgs[org] += 1
            #End of OG: Computa mean copye number
            Totcopy = 0
            for org in seen_orgs:
                number = seen_orgs[org]
                Totcopy = Totcopy + number
            orgnumber = len(seen_orgs)
            meanCopy = int(Totcopy)/int(orgnumber)
            if (meanCopy <= float(maxcopy)):
                tempOGlist_of.append(OG)
        #Replace list of OG
        totOGcount = len(OG_list)
        OG_list = tempOGlist_of
        keepOGcount = len(OG_list)
        lostOG = int(totOGcount) - int(keepOGcount)
        #Print stat
        string1= 'cFilter activated: Conserve only OG with ' + str(maxcopy) + ' in mean of sequence number of diffreent org.'
        out_file.write('#'+ str(string1) + "\n")
        string2= str(keepOGcount) + ' conserved OGs and ' + str(lostOG) + ' lost OGs.'
        out_file.write('#' + str(string2) + "\n")
        print(string1)
        print(string2)
    
    #Homogeneity filter (only with anvio)
    if (hfilter == 'yes'):  
        tempOGlist_of = []
        Geo_index = {}
        Fun_index = {}
        #open hindex file
        indexfile = open(hindex)
        for line in indexfile:
            index_record = line.replace("\n", "")
            split_list = index_record.split("\t")
            OG = split_list[0]
            functional_homogeneity_index = split_list[2]
            geometric_homogeneity_index = split_list[3]
            Geo_index[OG]= geometric_homogeneity_index
            Fun_index[OG]= functional_homogeneity_index
        #Read OG and filter based in index
        for OG in OG_list:
            OG = OG.replace("_OG.fasta", "")
            Geo = Geo_index[OG]
            Fun = Fun_index[OG]
            if ((Geo >= maxgeometricindex) and (Fun >= maxfunctionalindex)):
                string = str(OG) + '_OG.fasta'
                tempOGlist_of.append(string)
        #replace OG list
        totOGcount = len(OG_list)
        OG_list = tempOGlist_of
        keepOGcount = len(OG_list)
        lostOG = int(totOGcount) - int(keepOGcount)
        #Print stat
        string1= 'hFilter activated: Conserve only OG with ' + str(maxgeometricindex) + ' in geometric index and ' + str(maxfunctionalindex) + ' in functional index.' 
        string2= str(keepOGcount) + ' conserved OGs and ' + str(lostOG) + ' lost OGs.'
        out_file.write('#' + str(string1) + "\n")
        out_file.write('#' + str(string2) + "\n")
        print(string1)
        print(string2)
    
    #Read COG table and store COG 
    COG20_PATHWAY = {}
    COG20_FUNCTION = {}
    COG20_CATEGORY = {}
    #open cog table
    cogfile = open(cogtable)
    for line in cogfile:
        cog_record = line.replace("\n", "")
        split_list = cog_record.split("\t")
        OG = split_list[1]
        Pathway = split_list[13]
        Function = split_list[15]
        Category = split_list[17]
        COG20_PATHWAY[OG] = Pathway
        COG20_FUNCTION[OG] = Function
        COG20_CATEGORY[OG] = Category


    #print list in file
    #out_file = open('filtered-OG.list', "w")
    #print config
    out_file.write('#file:' + str(main_file) + "\n")
    out_file.write('#mode:' + str(mode) + "\n")
    out_file.write('#pfilter:' + str(pfilter) + "\n")
    out_file.write('#fraction:' + str(fraction) + "\n")
    out_file.write('#unwanted:' + str(unwanted) + "\n")
    out_file.write('#cfilter:' + str(cfilter) + "\n")
    out_file.write('#maxcopy:' + str(maxcopy) + "\n")
    out_file.write('#hfilter:' + str(hfilter) + "\n")
    out_file.write('#hindex:' + str(hindex) + "\n")
    out_file.write('#maxfunctionalindex:' + str(maxfunctionalindex) + "\n")
    out_file.write('#maxgeometricindex:' + str(maxgeometricindex) + "\n")
    out_file.write('OG-name' + "\t" + 'COG20_PATHWAY' + "\t" + 'COG20_FUNCTION' + "\t" + 'COG20_CATEGORY' + "\n")

    for OG in OG_list:
        OG = OG.replace("_OG.fasta", "")
        Pathway = COG20_PATHWAY[OG]
        Function = COG20_FUNCTION[OG]
        Category = COG20_CATEGORY[OG]
        out_file.write(str(OG) + "\t" + '|' + str(Pathway) + '|' + "\t" + '|' + str(Function) + '|' + "\t" + '|' + str(Category) + '|' + "\n")



if __name__ == '__main__':
    main()

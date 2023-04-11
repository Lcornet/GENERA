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
@click.option('--mode', default='no', help='rnammer or checkm mode')
@click.option('--ssu', default='no', help='delete genomes that without ssu : yes or no')
@click.option('--taxa', default='no', help='Choice between four taxa levels: phylum, class, order, family')
@click.option('--refgroup', default='no', help='Name of the refence group')
@click.option('--shrinkvalue', default='0.5', help='Three shrink cutoff value')

def main(main_file, mode, ssu, taxa, refgroup, shrinkvalue):
		
	#read fasta file

    if (mode == 'rnammer'):
        #Read ssu file and register GCF
        GCF_of = {}
        infile_ssu = open(main_file)
        rnammer_out = open("genome-with-ssu.list", "w")
        SSU_out = open("all_16s-nodupe.fna", "w")
        check_print = 'disabled'

        for line in infile_ssu:
            ssu_record = line.replace("\n", "")
            if '>' in ssu_record:
                #By default activated check print
                check_print = 'activated'
                #In the deflines-get GCF
                split_list = ssu_record.split("|")
                GCF = split_list[0]
                GCF = GCF.replace(">rRNA_", "")
                #skip the rest and don't store/print if GCF already seen
                if GCF in GCF_of:
                    check_print = 'disabled'
                #store result
                if (check_print == 'activated'):
                    GCF_of[GCF] = 1
                    SSU_out.write('>' + GCF + "\n")
            else:
                #print sequence if check print activated
                if (check_print == 'activated'):
                    SSU_out.write(ssu_record + "\n")

        #Print results 
        for GCF in GCF_of:
            rnammer_out.write(GCF + "\n")

    if (mode == 'barnap'):
        #Read ssu file and register GCF
        GCF_of = {}
        infile_ssu = open(main_file)
        rnammer_out = open("genome-with-ssu.list", "w")
        SSU_out = open("all_16s-nodupe.fna", "w")
        check_print = 'disabled'

        for line in infile_ssu:
            ssu_record = line.replace("\n", "")
            if ('>' in ssu_record) :
                #By default activated check print
                check_print = 'disabled'
                if ('16S' in ssu_record):
                    check_print = 'activated'
                #In the deflines-get GCF
                split_list = ssu_record.split("|")
                GCF = split_list[0]
                GCF = GCF.replace(">16S_rRNA::", "")
                #skip the rest and don't store/print if GCF already seen
                if GCF in GCF_of:
                    check_print = 'disabled'
                #store result
                if (check_print == 'activated'):
                    GCF_of[GCF] = 1
                    SSU_out.write('>' + GCF + "\n")
            else:
                #print sequence if check print activated
                if (check_print == 'activated'):
                    SSU_out.write(ssu_record + "\n")


        #Print results 
        for GCF in GCF_of:
            rnammer_out.write(GCF + "\n")
    
    elif (mode == 'checkm'):
        #define
        file_lists = []

        #outfile
        checkm_out = open("reliable-genomes.list", "w")

        #Open checkm result
        infile = open(main_file)

        for line in infile:
            checkm_record = line.replace("\n", "")
        
            #header
            if '#' in checkm_record:
                continue

            #split
            split_list = checkm_record.split(",")
            GCF = split_list[0]
            GCF = GCF.replace("-abbr", "")
            Compl = float(split_list[1])
            Conta = float(split_list[2])

            #criteria of filtration
            if ((Compl > 90) and (Conta < 5)):
                file_lists.append(GCF)
        
        #print result
        if (ssu == 'yes'):
            #delete genomes with ssu
            infile_rna = open("genome-with-ssu.list")
            for GCF in infile_rna:
                GCF_record = GCF.replace("\n", "")
                print(GCF_record)
                if GCF_record in file_lists:
                    checkm_out.write(GCF_record + "\n")

        elif (ssu == 'no'):
            #don't care of ssu
            for GCF in file_lists:
                checkm_out.write(GCF + "\n")

    elif (mode=='forty'):
        #define
        file_lists = []

        #outfile
        forty_out = open("enrich.list", "w")

        #Open checkm result
        infile = open(main_file)

        #read count file
        for line in infile:
            forty_record = line.replace("\n", "")
            split_list = forty_record.split(":")
            name = split_list[0]
            name = name.replace(".fasta", "")
            count = int(split_list[1])
            if (count >= 2):
                file_lists.append(name)

        #print results
        for name in file_lists:
            forty_out.write(name + "\n")
    
    elif (mode == 'sum'):
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
                    #print(sum_record)
                    split_list = sum_record.split("\t")
                    link = split_list[19]
                    if (link != 'na'):
                        sum_out.write(sum_record + "\n")
    
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
            
            #print GC list
            fetch_out = open("GCA.refgroup.uniq", "w")
            for GC in GCAlist_of:
                fetch_out.write(str(GC) + "\n")

    elif (mode == 'ConstrainTreeRaxml'):
        
        #open list of GC ids
        GC_infile = open('GC.list')
        GClist_of = {}
        for line in GC_infile:
            GC = line.replace("\n", "")
            GClist_of[GC] = 1

        #Open refence 16s from rnammer
        infile = open(main_file)

        #Open out file
        GC_out = open("all_16s-nodupe-list.fasta", "w")
        
        tag = 0
        for line in infile:
            GC_record = line.replace("\n", "")
            #activate tag in GC in list
            if ('>' in GC_record):
                #reset tag at each new id
                tag = 0
                #check
                GC = GC_record.replace(">", "")
                if (GC in GClist_of):
                    tag += 1
            #print if tag activated
            if (tag == 1):
                GC_out.write(str(GC_record) + "\n")
    
    elif (mode =='shrink'):

        #produce GCA list of out group
        ref_GCA = {}
        reffile = open('Out-GCA.list')
        for line in reffile:
            ref_record = line.replace("\n", "")
            ref_GCA[ref_record]= 1

        #open file from treeshrinl
        infile = open(main_file)

        #print idl
        idl_out = open("Constained-tree.idl", "w")
        idl_out.write('#' + "\n")
        idl_out.write('#' + "\n")

        #pass into treeshrink file and produce idl
        threshold = shrinkvalue
        for line in infile:
            tag = 1
            idl_record = line.replace("\n", "")
            #check if line is not about outgroup
            for ref in ref_GCA:
                if ref in idl_record:
                    tag = 0
            if 'Signature' in idl_record:
                continue
            else:
                split_list = idl_record.split(" ")
                treesh_value = float(split_list[3])
                org = split_list[2]
                if ((treesh_value >= float(threshold)) and (tag == 1)):
                    idl_out.write(str(org) + "\n")




            

if __name__ == '__main__':
    main()
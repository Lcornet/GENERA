#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Companion"""
__author__ = "Luc Cornet"
__copyright__ = "Copyright 2021, Sciensano"
__version__ = "1.0.0"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"

from unicodedata import category
import click
import glob
import random
from Bio.Seq import Seq
from collections import defaultdict

@click.command()
###ARGUMENT#####
#main file - depending of the mode can be checkm result or rnammer output
@click.argument('main_file', type=click.Path(exists=True,readable=True))
#taxid
@click.option('--mode', default='no', help='mode')
@click.option('--mask_slave', default='no', help='mode')
@click.option('--level', default='no', help='mode')
@click.option('--number', default='100', help='mode')
@click.option('--corr_file', default='no')
@click.option('--redundancy', default='20')
@click.option('--replacement', default='20')
@click.option('--single', default='20')
@click.option('--merge', default='no', help='merge')
@click.option('--redundancyhgt', default='20')
@click.option('--replacementhgt', default='20')
@click.option('--singlehgt', default='20')

def main(main_file, mode, mask_slave, level, number, corr_file, redundancy, replacement, single, merge, redundancyhgt, replacementhgt, singlehgt):

    #read taxonomy file
    GC_taxoDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
    GC_list = {}
    if (mode == 'list'):
        #Open main file
        infile = open(main_file)
        for line in infile:
            record = line.replace("\n", "")
            split_list = record.split("\t")
            GC = split_list[0]
            name = split_list[1]
            lineage = split_list[3]
            split_list = lineage.split(";")
            phylum_r = split_list[0]
            class_r = split_list[1] 
            order_r = split_list[2] 
            family_r = split_list[3] 
            genus_r = split_list[4]
            species_r = split_list[5]
            #Load in nested dict
            GC_taxoDD[GC]['name'] = name
            GC_taxoDD[GC]['phylum'] = phylum_r
            GC_taxoDD[GC]['class'] = class_r
            GC_taxoDD[GC]['order'] = order_r
            GC_taxoDD[GC]['family'] = family_r
            GC_taxoDD[GC]['genus'] = genus_r
            GC_taxoDD[GC]['species'] = species_r
            GC_list[GC] = 1
    
        #Create the corresponding list
        list_out = open('chimeric.list', "w")
        corr_out = open('chimeric.corr', "w")
        out = 'chimeric.idl'
        idl_out = open(out, "w")

        rank_used = {}

        if (level == 'interphylum'):
            used_GC = {}
            hits = 0
            #Shuffle the dict
            list_ofs = list(GC_list.items())
            random.shuffle(list_ofs)
            GC_shuff = dict(list_ofs)
            for master in GC_shuff:
                if (hits == int(number)):
                    break
                #Get master rank
                masterRank1 = GC_taxoDD[master]['phylum'] 
                if (masterRank1 == 'undef'):
                    continue
                #Second loop: slave
                for slave in GC_shuff:
                    if ((slave == master) or (slave in used_GC)):
                        continue
                    else:
                        #Get slave rank
                        slaveRank1 = GC_taxoDD[slave]['phylum']
                        if (slaveRank1 == 'undef'):
                            continue
                        #The two GC should of the same taxonomic rank1
                        #The two GC should not be of the same taxonomic rank2
                        if (masterRank1 != slaveRank1):
                            #print(masterRank2 + " - " + slaveRank2)
                            #corresponding match
                            used_GC[master]=1
                            if (mask_slave == 'yes'):
                                used_GC[slave]=1
                            hits +=1
                            if (masterRank1 not in rank_used):
                                rank_used[masterRank1] = 1
                                idl_out.write(str(masterRank1) + "\n")
                            if (slaveRank1 not in rank_used):
                                rank_used[slaveRank1] = 1
                                idl_out.write(str(slaveRank1) + "\n")                   
                            corr_out.write(str(master) + "\t" + str(masterRank1) + "\t" + str(slave) + "\t" + str(slaveRank1) +  "\n")
                            list_out.write(str(master) + "\n")
                            list_out.write(str(slave) + "\n")
                            break

        if (level == 'phylum'):
            used_GC = {}
            hits = 0
            #Shuffle the dict
            list_ofs = list(GC_list.items())
            random.shuffle(list_ofs)
            GC_shuff = dict(list_ofs)
            for master in GC_shuff:
                if (hits == int(number)):
                    break
                #Get master rank
                masterRank1 = GC_taxoDD[master]['phylum'] 
                masterRank2 = GC_taxoDD[master]['class']
                if ((masterRank1 == 'undef') or (masterRank2 == 'undef')):
                    continue
                #Second loop: slave
                for slave in GC_shuff:
                    if ((slave == master) or (slave in used_GC)):
                        continue
                    else:
                        #Get slave rank
                        slaveRank1 = GC_taxoDD[slave]['phylum']
                        slaveRank2 = GC_taxoDD[slave]['class']
                        if ((slaveRank1 == 'undef') or (slaveRank2 == 'undef')):
                            continue
                        #The two GC should of the same taxonomic rank1
                        #The two GC should not be of the same taxonomic rank2
                        if ((masterRank1 == slaveRank1) and (masterRank2 != slaveRank2)):
                            #print(masterRank2 + " - " + slaveRank2)
                            #corresponding match
                            used_GC[master]=1
                            if (mask_slave == 'yes'):
                                used_GC[slave]=1
                            hits +=1
                            if (masterRank2 not in rank_used):
                                rank_used[masterRank2] = 1
                                idl_out.write(str(masterRank2) + "\n")
                            if (slaveRank2 not in rank_used):
                                rank_used[slaveRank2] = 1
                                idl_out.write(str(slaveRank2) + "\n")                   
                            corr_out.write(str(master) + "\t" + str(masterRank1) + '-' + str(masterRank2) + "\t" + str(slave) + "\t" + str(slaveRank1) + '-' + str(slaveRank2) + "\n")
                            list_out.write(str(master) + "\n")
                            list_out.write(str(slave) + "\n")
                            break
        
        if (level == 'class'):
            used_GC = {}
            hits = 0
            #Shuffle the dict
            list_ofs = list(GC_list.items())
            random.shuffle(list_ofs)
            GC_shuff = dict(list_ofs)
            for master in GC_shuff:
                if (hits == int(number)):
                    break
                #Get master rank
                masterRank1 = GC_taxoDD[master]['class'] 
                masterRank2 = GC_taxoDD[master]['order']
                if ((masterRank1 == 'undef') or (masterRank2 == 'undef')):
                    continue
                #Second loop: slave
                for slave in GC_shuff:
                    if ((slave == master) or (slave in used_GC)):
                        continue
                    else:
                        #Get slave rank
                        slaveRank1 = GC_taxoDD[slave]['class']
                        slaveRank2 = GC_taxoDD[slave]['order']
                        if ((slaveRank1 == 'undef') or (slaveRank2 == 'undef')):
                            continue
                        #The two GC should of the same taxonomic rank1
                        #The two GC should not be of the same taxonomic rank2
                        if ((masterRank1 == slaveRank1) and (masterRank2 != slaveRank2)):
                            #print(masterRank2 + " - " + slaveRank2)
                            #corresponding match
                            used_GC[master]=1
                            if (mask_slave == 'yes'):
                                used_GC[slave]=1
                            hits +=1
                            if (masterRank2 not in rank_used):
                                rank_used[masterRank2] = 1
                                idl_out.write(str(masterRank2) + "\n")
                            if (slaveRank2 not in rank_used):
                                rank_used[slaveRank2] = 1
                                idl_out.write(str(slaveRank2) + "\n")                   
                            corr_out.write(str(master) + "\t" + str(masterRank1) + '-' + str(masterRank2) + "\t" + str(slave) + "\t" + str(slaveRank1) + '-' + str(slaveRank2) + "\n")
                            list_out.write(str(master) + "\n")
                            list_out.write(str(slave) + "\n")
                            break      

        if (level == 'order'):
            used_GC = {}
            hits = 0
            #Shuffle the dict
            list_ofs = list(GC_list.items())
            random.shuffle(list_ofs)
            GC_shuff = dict(list_ofs)
            for master in GC_shuff:
                if (hits == int(number)):
                    break
                #Get master rank
                masterRank1 = GC_taxoDD[master]['order'] 
                masterRank2 = GC_taxoDD[master]['family']
                if ((masterRank1 == 'undef') or (masterRank2 == 'undef')):
                    continue
                #Second loop: slave
                for slave in GC_shuff:
                    if ((slave == master) or (slave in used_GC)):
                        continue
                    else:
                        #Get slave rank
                        slaveRank1 = GC_taxoDD[slave]['order']
                        slaveRank2 = GC_taxoDD[slave]['family']
                        if ((slaveRank1 == 'undef') or (slaveRank2 == 'undef')):
                            continue
                        #The two GC should of the same taxonomic rank1
                        #The two GC should not be of the same taxonomic rank2
                        if ((masterRank1 == slaveRank1) and (masterRank2 != slaveRank2)):
                            #print(masterRank2 + " - " + slaveRank2)
                            #corresponding match
                            used_GC[master]=1
                            if (mask_slave == 'yes'):
                                used_GC[slave]=1
                            hits +=1
                            if (masterRank2 not in rank_used):
                                rank_used[masterRank2] = 1
                                idl_out.write(str(masterRank2) + "\n")
                            if (slaveRank2 not in rank_used):
                                rank_used[slaveRank2] = 1
                                idl_out.write(str(slaveRank2) + "\n")  
                            corr_out.write(str(master) + "\t" + str(masterRank1) + '-' + str(masterRank2) + "\t" + str(slave) + "\t" + str(slaveRank1) + '-' + str(slaveRank2) + "\n")
                            list_out.write(str(master) + "\n")
                            list_out.write(str(slave) + "\n")
                            break

        if (level == 'family'):
            used_GC = {}
            hits = 0
            #Shuffle the dict
            list_ofs = list(GC_list.items())
            random.shuffle(list_ofs)
            GC_shuff = dict(list_ofs)
            for master in GC_shuff:
                if (hits == int(number)):
                    break
                #Get master rank
                masterRank1 = GC_taxoDD[master]['family'] 
                masterRank2 = GC_taxoDD[master]['genus']
                if ((masterRank1 == 'undef') or (masterRank2 == 'undef')):
                    continue
                #Second loop: slave
                for slave in GC_shuff:
                    if ((slave == master) or (slave in used_GC)):
                        continue
                    else:
                        #Get slave rank
                        slaveRank1 = GC_taxoDD[slave]['family']
                        slaveRank2 = GC_taxoDD[slave]['genus']
                        if ((slaveRank1 == 'undef') or (slaveRank2 == 'undef')):
                            continue
                        #The two GC should of the same taxonomic rank1
                        #The two GC should not be of the same taxonomic rank2
                        if ((masterRank1 == slaveRank1) and (masterRank2 != slaveRank2)):
                            #print(masterRank2 + " - " + slaveRank2)
                            #corresponding match
                            used_GC[master]=1
                            if (mask_slave == 'yes'):
                                used_GC[slave]=1
                            hits +=1
                            if (masterRank2 not in rank_used):
                                rank_used[masterRank2] = 1
                                idl_out.write(str(masterRank2) + "\n")
                            if (slaveRank2 not in rank_used):
                                rank_used[slaveRank2] = 1
                                idl_out.write(str(slaveRank2) + "\n")  
                            corr_out.write(str(master) + "\t" + str(masterRank1) + '-' + str(masterRank2) + "\t" + str(slave) + "\t" + str(slaveRank1) + '-' + str(slaveRank2) + "\n")
                            list_out.write(str(master) + "\n")
                            list_out.write(str(slave) + "\n")
                            break
    

        if (level == 'genus'):
            used_GC = {}
            hits = 0
            #Shuffle the dict
            list_ofs = list(GC_list.items())
            random.shuffle(list_ofs)
            GC_shuff = dict(list_ofs)
            for master in GC_shuff:
                if (hits == int(number)):
                    break
                #Get master rank
                masterRank1 = GC_taxoDD[master]['genus'] 
                masterRank2 = GC_taxoDD[master]['species']
                if ((masterRank1 == 'undef') or (masterRank2 == 'undef')):
                    continue
                #Second loop: slave
                for slave in GC_shuff:
                    if ((slave == master) or (slave in used_GC)):
                        continue
                    else:
                        #Get slave rank
                        slaveRank1 = GC_taxoDD[slave]['genus']
                        slaveRank2 = GC_taxoDD[slave]['species']
                        if ((slaveRank1 == 'undef') or (slaveRank2 == 'undef')):
                            continue
                        #The two GC should of the same taxonomic rank1
                        #The two GC should not be of the same taxonomic rank2
                        if ((masterRank1 == slaveRank1) and (masterRank2 != slaveRank2)):
                            #print(masterRank2 + " - " + slaveRank2)
                            #corresponding match
                            used_GC[master]=1
                            if (mask_slave == 'yes'):
                                used_GC[slave]=1
                            hits +=1
                            if (masterRank2 not in rank_used):
                                rank_used[masterRank2] = 1
                                idl_out.write(str(masterRank2) + "\n")
                            if (slaveRank2 not in rank_used):
                                rank_used[slaveRank2] = 1
                                idl_out.write(str(slaveRank2) + "\n")   
                            corr_out.write(str(master) + "\t" + str(masterRank1) + '-' + str(masterRank2) + "\t" + str(slave) + "\t" + str(slaveRank1) + '-' + str(slaveRank2) + "\n")
                            list_out.write(str(master) + "\n")
                            list_out.write(str(slave) + "\n")
                            break


        if (level == 'species'):
            used_GC = {}
            hits = 0
            #Shuffle the dict
            list_ofs = list(GC_list.items())
            random.shuffle(list_ofs)
            GC_shuff = dict(list_ofs)
            for master in GC_shuff:
                if (hits == int(number)):
                    break
                #Get master rank
                masterRank1 = GC_taxoDD[master]['species'] 
                if (masterRank1 == 'undef'):
                    continue
                #Second loop: slave
                for slave in GC_shuff:
                    if ((slave == master) or (slave in used_GC)):
                        continue
                    else:
                        #Get slave rank
                        slaveRank1 = GC_taxoDD[slave]['species']
                        if (slaveRank1 == 'undef'):
                            continue
                        #The two GC should of the same taxonomic rank1
                        if (masterRank1 == slaveRank1):
                            #print(masterRank1 + " - " + slaveRank1)
                            #corresponding match
                            used_GC[master]=1
                            if (mask_slave == 'yes'):
                                used_GC[slave]=1
                            hits +=1
                            if (masterRank1 not in rank_used):
                                rank_used[masterRank1] = 1
                                idl_out.write(str(masterRank1) + "\n")
                            if (slaveRank1 not in rank_used):
                                rank_used[slaveRank1] = 1
                                idl_out.write(str(slaveRank1) + "\n") 
                            corr_out.write(str(master) + "\t" + str(masterRank1) + "\t" + str(slave) + "\t" + str(slaveRank1) +  "\n")
                            list_out.write(str(master) + "\n")
                            list_out.write(str(slave) + "\n")
                            break

    if (mode == 'chimeric'):

        chim_log = open('chimeric.log', "w")

        genome_list = {}
        #Open List of genomes involved
        infile = open(main_file)
        for line in infile:
            record = line.replace("\n", "")
            genome_list[record] = 1

        #Load sequences of genomes
        GC_fnaDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        GC_cooDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        Genome_num = {}
        for genome in genome_list:
            genes_file = str(genome) + ".genes.fna"
            defline = 'NA'
            deflineNum = 0
            sequence = None
            start = 'NA'
            end = 'NA'
            ID = 'NA'
            strand = 'NA'
            #file reading
            genes = open(genes_file)
            for line in genes:
                record = line.replace("\n", "")
                if '>' in line:
                    #fasta reading mode
                    if (defline != 'NA'):
                        #load the genes and reset
                        GC_fnaDD[genome][defline] = sequence
                        GC_cooDD[genome][defline]['start'] = start
                        GC_cooDD[genome][defline]['end'] = end
                        GC_cooDD[genome][defline]['ID'] = ID
                        GC_cooDD[genome][defline]['strand'] = strand
                        defline = 'NA'
                        sequence = None
                        start = 'NA'
                        end = 'NA'
                        ID = 'NA'
                        strand = 'NA'
                    split_list = record.split(" ")
                    defline = split_list[0]
                    start = split_list[2]
                    end = split_list[4]
                    strand = split_list[6]
                    ID = split_list[8]
                    split_list = ID.split(";")
                    ID = split_list[0]
                    split_list = ID.split("_")
                    ID = split_list[1]
                    defline = defline.replace(">", "")
                    split_list = defline.split("|")
                    defline = split_list[1]
                    deflineNum += 1
                else:
                    #cat sequences
                    if (sequence == None):
                        sequence = str(record)
                    else:
                        sequence += str(record)
            #final sequence load
            GC_fnaDD[genome][defline] = sequence
            GC_cooDD[genome][defline]['start'] = start
            GC_cooDD[genome][defline]['end'] = end
            GC_cooDD[genome][defline]['ID'] = ID
            GC_cooDD[genome][defline]['strand'] = strand
            Genome_num[genome] = deflineNum

        #load full sequences from contigs
        GC_genoDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        GC_corDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        for genome in genome_list:
            genome_file = str(genome) + ".fna"
            defline = 'NA'
            sequence = None
            #files reading
            geno = open(genome_file)
            for line in geno:
                record = line.replace("\n", "")
                if '>' in line:
                    #fasta reading mode
                    if (defline != 'NA'):
                        #load the genes and reset
                        GC_genoDD[genome][defline] = sequence
                        GC_corDD[genome][defline] = 0
                        defline = 'NA'
                        sequence = None
                    split_list = record.split(" ")
                    defline = split_list[0]
                    defline = defline.replace(">", "")
                    split_list = defline.split("|")
                    defline = split_list[1]
                else:
                    #cat sequences
                    if (sequence == None):
                        sequence = str(record)
                    else:
                        sequence += str(record)
            #final sequence load
            GC_genoDD[genome][defline] = sequence
            GC_corDD[genome][defline] = 0

        #print(GC_corDD)

        #Open OG files and exclude OG with multi copy genes
        OG_list = glob.glob("OG/OG*.fa")
        OG_excluded = {}
        OG_deflinesDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))

        for OG in OG_list:
            OGfile = open(OG)
            genome_seen = {}
            for line in OGfile:
                record = line.replace("\n", "")
                if '>' in line:
                    defline = record.replace(">", "")
                    split_list = defline.split("|")
                    genome = split_list[0]
                    if (genome in genome_seen):
                        OG_excluded[OG] = 1
                    else:
                        genome_seen[genome] = 1
                    gene = split_list[1]
                    #load deflines
                    OG_deflinesDD[OG][genome] = gene

        #Make replacement
        #Load corr.list
        corrfile = open(corr_file)
        corr_list = {}
        taxo_list = {}
        for line in corrfile:
            record = line.replace("\n", "")
            split_list = record.split("\t")
            genomeM = split_list[0]
            taxoM = split_list[1]
            genomeS = split_list[2]
            taxoS = split_list[3]
            corr_list[genomeM] = genomeS
            taxo_list[genomeM] = taxoM
            taxo_list[genomeS] = taxoS

        #make chimeric
        genumber = 0
        for master in corr_list:
            #Number of defline
            Msize = Genome_num[master]
            Replacement_size = int((float(Msize) / 100) * float(replacement))
            Duplication_size = int((float(Msize) / 100) * float(redundancy))
            Single_size = int((float(Msize) / 100) * float(single))
            ReplacementHGT_size = int((float(Msize) / 100) * float(replacementhgt))
            DuplicationHGT_size = int((float(Msize) / 100) * float(redundancyhgt))
            SingleHGT_size = int((float(Msize) / 100) * float(singlehgt))
            DupREp_size = Replacement_size + Duplication_size
            DupRepHGT_size = ReplacementHGT_size + DuplicationHGT_size
            DupREpALL_size = DupREp_size + DupRepHGT_size
            SingleALL_size = Single_size + SingleHGT_size
            #Change_size = Replacement_size + Duplication_size + Single_size
            #slave 
            slave = corr_list[master]
            
            #Check if redudancy and replacement replacement are possible
            #pass in OGs and count number of DupREp presence in single copy genes
            common_OG = {}
            single_OG = {}
            for OG in OG_list:
                #skip non single copy OG
                if (OG in OG_excluded):
                    continue
                #def
                masterSeen = 0
                slaveSeen = 0
                #denest and count
                genome_of = OG_deflinesDD[OG]
                for genome in genome_of:
                    if (genome == master):
                        masterSeen += 1
                    elif (genome == slave):
                        slaveSeen += 1
                #Do the OG pass ? for redundancy and replacement
                if (masterSeen == 1) and (slaveSeen == 1):
                    common_OG[OG] = 1
                #Do the OG pass ? for single
                if (masterSeen == 0) and (slaveSeen == 1):
                    single_OG[OG] = 1     

            #Delete overlapping genes
            OGstart_orderedDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
            OGend_orderedDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
            for OG in common_OG:
                genome_of = OG_deflinesDD[OG]
                for genome in genome_of:
                    if (genome == master):
                        masterGene = genome_of[master]
                        master_ID = '_' + str(GC_cooDD[master][masterGene]['ID'])
                        master_contig = masterGene.replace(master_ID, "")
                        master_start = int(GC_cooDD[master][masterGene]['start'])
                        master_end = int(GC_cooDD[master][masterGene]['end'])
                        OGstart_orderedDD[master_contig][OG]= master_start
                        OGend_orderedDD[master_contig][OG]= master_end
            #Delete gene if coordinate overlap previous one
            previous_end = 'NA'
            previous_OG = 'NA'
            for master_contig in OGstart_orderedDD:
                contig_of = OGstart_orderedDD[master_contig]
                for OG in sorted(contig_of, key=contig_of.get, reverse=False):
                    start = OGstart_orderedDD[master_contig][OG]
                    end = OGend_orderedDD[master_contig][OG]
                    if (previous_end == 'NA'):
                        previous_end = end
                        previous_OG = OG
                    else:
                        if (int(start) < int(previous_end)):
                            del common_OG[OG]
                            if (previous_OG in common_OG):
                                del common_OG[previous_OG]
                            previous_end = end
                            previous_OG = OG
                        else:
                            previous_end = end
                            previous_OG = OG


            #Final check
            commonOG = len(common_OG)
            singelOG = len(single_OG)
            availableOG = commonOG + singelOG
            taxoM = taxo_list[master]
            taxoS = taxo_list[slave]

            #if (Change_size > availableOG):
            if ((float(DupREpALL_size) > float(commonOG)) or (float(SingleALL_size) > float(singelOG))):
                #cannot create
                chim_log.write(str(master) + "\t" + str(taxoM) + "\t" + str(slave) + "\t" + str(taxoS)
                + "\t" +  "The number of available proteins in Sub genome are not sufficient: "  
                + str(commonOG) + " proteins while " +  str(DupREpALL_size) + " are needed for redundancy and/or replacement - " 
                + str(singelOG) + " proteins while " +  str(SingleALL_size) + " are needed for singel." + "\n")
            else:
                #Create a hash to order the OGs according to their starting position on master genome
                commonOG_orderedDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
                #only load commonOG in, singleOG not integrated within master_ctg as replacement with commonOG
                #SingleOG treated separately
                for OG in common_OG:
                    genome_of = OG_deflinesDD[OG]
                    for genome in genome_of:
                        if (genome == master):
                            masterGene = genome_of[master]
                            master_ID = '_' + str(GC_cooDD[master][masterGene]['ID'])
                            master_contig = masterGene.replace(master_ID, "")
                            master_start = int(GC_cooDD[master][masterGene]['start'])
                            commonOG_orderedDD[master_contig][OG] = master_start 

                #Create chimeric genome
                genumber += 1
                #chim_log.write(str(master) + "\t" + str(slave) + "\t" + 'creation of chimeric genome: ID=' + " " + str(genumber) + "\n")

                genome_name = 'chimeric-' + str(genumber) + ".ali"
                gen_out = open(genome_name, "w")
                genome_log = 'chimeric-' + str(genumber) + "-log.ali"
                gen_log = open(genome_log, "w")
                HGTcoo_log = 'chimeric-' + str(genumber) + "-HGT.coo"
                hgtcoo_log = open(HGTcoo_log, "w")

                #Open file for HGTsim
                #Open two file for HGIsim, one with name of genome + genes to be transferred (Distribution), and on with sequence of genes to be transferred
                #HGTsim distribution
                HGTsimdistibution = 'chimeric-' + str(genumber) + '-distribution.txt'
                HGTsim_distribution = open(HGTsimdistibution, "w")
                HGTsim_distribution.write('chimeric-' + str(genumber))
                #HGTsim genes
                HGTsimgenes = 'chimeric-' + str(genumber) + '-genes.ali'
                HGTsim_genes = open(HGTsimgenes, "w")


                #Pass in master genome and process to replacement
                #First make replacement from OGs
                gene_done = {}
                #Define
                replacementyNumber = 0
                duplicationNumber = 0
                contam_size = 0
                chim_size = 0
                replacementyNumberHGT = 0
                duplicationNumberHGT = 0
                HGT_size = 0             
                NEW = 1
                #Pass into common OG for Duplication and replacement
                duplication_ctg = None
                #for OG in common_OG:
                for master_contig in commonOG_orderedDD:
                    contig_of = commonOG_orderedDD[master_contig]
                    for OG in sorted(contig_of, key=contig_of.get, reverse=False):
                        #Denest
                        genome_of = OG_deflinesDD[OG]
                        masterGene = 'na'
                        slaveGene = 'na'
                        for genome in genome_of:
                            if (genome == master):
                                masterGene = genome_of[master]
                            elif (genome == slave):
                                slaveGene = genome_of[slave]
                        #Duplication first
                        if (duplicationNumber < Duplication_size):
                            #Print DupREp
                            #master
                            #For the genome: Nothing to do
                            #For the log
                            gene_done[masterGene] = 1
                            master_ID = '_' + str(GC_cooDD[master][masterGene]['ID'])
                            master_contig = masterGene.replace(master_ID, "")
                            master_ctg = GC_genoDD[master][master_contig]
                            master_cor = GC_corDD[master][master_contig]
                            master_start = int(GC_cooDD[master][masterGene]['start']) - 1 + int(master_cor)
                            master_end = int(GC_cooDD[master][masterGene]['end']) + int(master_cor)
                            #slave
                            slave_ID = '_' + str(GC_cooDD[slave][slaveGene]['ID'])
                            slave_contig = slaveGene.replace(slave_ID, "")
                            slave_ctg = GC_genoDD[slave][slave_contig]
                            slave_start = int(GC_cooDD[slave][slaveGene]['start']) - 1
                            slave_end = int(GC_cooDD[slave][slaveGene]['end']) 
                            slave_seq = slave_ctg[slave_start:slave_end]
                            #For the genome
                            if (duplication_ctg == None):
                                duplication_ctg = 'NNNNN' + str(slave_seq)
                            else:
                                duplication_ctg = str(duplication_ctg) + 'NNNNN' + str(slave_seq)
                            #NEW += 1
                            contam_size += len(slave_seq)
                            #for the log
                            gen_log.write(">" + str(slave) + '|' + str(slaveGene) + '|' + 'Duplicate: ' 
                            + str(masterGene) + ":" + str(master_start) + '-' + str(master_end) +  "\n")
                            gen_log.write(str(slave_seq) + "\n")                       
                            #Increment
                            duplicationNumber += 1
                        elif (replacementyNumber < Replacement_size):
                            #Replacement in second
                            #Print only slave
                            #master
                            gene_done[masterGene] = 1
                            master_ID = '_' + str(GC_cooDD[master][masterGene]['ID'])
                            master_contig = masterGene.replace(master_ID, "")
                            master_ctg = GC_genoDD[master][master_contig]
                            master_cor = GC_corDD[master][master_contig]
                            master_start = int(GC_cooDD[master][masterGene]['start']) - 1 + int(master_cor)
                            master_end = int(GC_cooDD[master][masterGene]['end']) + int(master_cor)
                            master_seq = master_ctg[master_start:master_end]
                            #slave
                            #For the genome
                            slave_ID = '_' + str(GC_cooDD[slave][slaveGene]['ID'])
                            slave_contig = slaveGene.replace(slave_ID, "")
                            slave_ctg = GC_genoDD[slave][slave_contig]
                            slave_start = int(GC_cooDD[slave][slaveGene]['start']) - 1
                            slave_end = int(GC_cooDD[slave][slaveGene]['end'])
                            slave_seq = slave_ctg[slave_start:slave_end]
                            #Compute length difference
                            mod_size = len(slave_seq) - len(master_seq)
                            mod_cor = int(master_cor) + int(mod_size)
                            GC_corDD[master][master_contig] = mod_cor
                            #Replacement of sequence in the hash
                            corrected_end = master_end + mod_cor
                            chim_seq = master_ctg[:master_start] + str(slave_seq) + master_ctg[master_end:]
                            GC_genoDD[master][master_contig] = chim_seq
                            #for the log
                            contam_size += len(slave_seq)
                            #For the log
                            gen_log.write(">" + str(slave) + '|' + str(slaveGene) + '|' + 'Replace: ' + str(masterGene) 
                            + ":" + str(master_start) + "-" + str(corrected_end) +  "\n")
                            gen_log.write(str(slave_seq) + "\n")                         
                            #Increment
                            replacementyNumber += 1
                        elif (duplicationNumberHGT < DuplicationHGT_size):
                            #Duplication of HGT in third
                            #master
                            #For the genome: Nothing to do
                            #For the log
                            gene_done[masterGene] = 1
                            master_ID = '_' + str(GC_cooDD[master][masterGene]['ID'])
                            master_contig = masterGene.replace(master_ID, "")
                            master_ctg = GC_genoDD[master][master_contig]
                            master_cor = GC_corDD[master][master_contig]
                            master_start = int(GC_cooDD[master][masterGene]['start']) - 1 + int(master_cor)
                            master_end = int(GC_cooDD[master][masterGene]['end']) + int(master_cor)
                            #slave
                            #For HGTsim, need only strand == 1
                            slave_ID = '_' + str(GC_cooDD[slave][slaveGene]['ID'])
                            slave_contig = slaveGene.replace(slave_ID, "")
                            slave_ctg = GC_genoDD[slave][slave_contig]
                            slave_start = int(GC_cooDD[slave][slaveGene]['start']) - 1
                            slave_end = int(GC_cooDD[slave][slaveGene]['end']) 
                            slave_seq = slave_ctg[slave_start:slave_end]
                            slave_strand = GC_cooDD[slave][slaveGene]['strand']
                            if (slave_strand == '-1'):
                                seq = Seq(slave_seq)
                                reverse = seq.reverse_complement()
                                slave_seq = reverse
                            #For the genome
                            #Print gene name in the list for HGTsim and gene seq in a separate file for HGTsim
                            HGTsim_distribution.write(',' + str(slaveGene))
                            HGTsim_genes.write('>' + str(slaveGene) + "\n")
                            HGTsim_genes.write(str(slave_seq) + "\n")
                            #NEW += 1
                            HGT_size += len(slave_seq)
                            #for the log
                            gen_log.write(">" + str(slave) + '|' + str(slaveGene) + '|' + 'DuplicateHGT-beforeMutations: ' 
                            + str(masterGene) + ":" + str(master_start) + '-' + str(master_end) +  "\n")
                            gen_log.write(str(slave_seq) + "\n")  
                            #For HGT coor 
                            hgtcoo_log.write('chimeric-' + str(genumber) + "\t" + str(master) + "\t" +  str(slaveGene) 
                            + "\t" + 'Duplication' + "\t" + str(master_contig) + "\t" + str(master_start) + "\t" + str(master_end) + "\n")                  
                            #Increment
                            duplicationNumberHGT += 1
                        elif (replacementyNumberHGT < ReplacementHGT_size):
                            #Replacement of HGT in fourth
                            #Print only slave
                            #master
                            gene_done[masterGene] = 1
                            master_ID = '_' + str(GC_cooDD[master][masterGene]['ID'])
                            master_contig = masterGene.replace(master_ID, "")
                            master_ctg = GC_genoDD[master][master_contig]
                            master_cor = GC_corDD[master][master_contig]
                            master_start = int(GC_cooDD[master][masterGene]['start']) - 1 + int(master_cor)
                            master_end = int(GC_cooDD[master][masterGene]['end']) + int(master_cor)
                            master_seq = master_ctg[master_start:master_end]
                            #slave
                            #For the genome
                            #Print gene name in the list for HGTsim and gene seq in a separate file for HGTsim
                            slave_ID = '_' + str(GC_cooDD[slave][slaveGene]['ID'])
                            slave_contig = slaveGene.replace(slave_ID, "")
                            slave_ctg = GC_genoDD[slave][slave_contig]
                            slave_start = int(GC_cooDD[slave][slaveGene]['start']) - 1
                            slave_end = int(GC_cooDD[slave][slaveGene]['end'])
                            slave_seq = slave_ctg[slave_start:slave_end]
                            slave_strand = GC_cooDD[slave][slaveGene]['strand']
                            if (slave_strand == '-1'):
                                seq = Seq(slave_seq)
                                reverse = seq.reverse_complement()
                                slave_seq = reverse                     
                            #Compute length difference
                            #Master gene is replaced by NNN..., not modification of length as it was for contamination replacement
                            mod_cor = int(master_cor) 
                            GC_corDD[master][master_contig] = mod_cor
                            #Replacement of sequence in the hash
                            master_length = len(master_seq)
                            Ntxt = ""
                            Nseq = Ntxt.ljust(master_length, "N")
                            chim_seq = master_ctg[:master_start] + str(Nseq) + master_ctg[master_end:]
                            GC_genoDD[master][master_contig] = chim_seq
                            #Outfile for HGTsim
                            HGTsim_distribution.write(',' + str(slaveGene))
                            HGTsim_genes.write('>' + str(slaveGene) + "\n")
                            HGTsim_genes.write(str(slave_seq) + "\n")
                            #for the log
                            HGT_size += len(slave_seq)
                            #For the log
                            gen_log.write(">" + str(slave) + '|' + str(slaveGene) + '|' + 'ReplaceHGT-beforeMutations: ' + str(masterGene) 
                            + ":" + str(master_start) + "-" + str(master_end) +  "\n")
                            gen_log.write(str(slave_seq) + "\n") 
                            #For HGT coor 
                            hgtcoo_log.write('chimeric-' + str(genumber) + "\t" + str(master) + "\t" +  str(slaveGene) 
                            + "\t" + 'Replace' + "\t" + str(master_contig) + "\t" + str(master_start) + "\t" + str(master_end) + "\n")                         
                            #Increment
                            replacementyNumberHGT += 1                       

                #Pass into singel OG, and create singel sequences
                singelNumber = 0
                singelNumberHGT = 0
                SIN = 0
                single_ctg = None
                for OG in single_OG:
                    #Denest
                    genome_of = OG_deflinesDD[OG]
                    slaveGene = 'na'
                    for genome in genome_of:
                        if (genome == slave):
                            slaveGene = genome_of[slave]
                    #Singel add
                    if (singelNumber < Single_size):
                        #Single contams in fifth
                        #For the genome: Nothing to do
                        #slave
                        slave_ID = '_' + str(GC_cooDD[slave][slaveGene]['ID'])
                        slave_contig = slaveGene.replace(slave_ID, "")
                        slave_ctg = GC_genoDD[slave][slave_contig]
                        slave_start = int(GC_cooDD[slave][slaveGene]['start']) - 1
                        slave_end = int(GC_cooDD[slave][slaveGene]['end'])
                        slave_seq = slave_ctg[slave_start:slave_end]
                        #For the genome
                        if (single_ctg == None):
                            single_ctg = 'NNNNN' + str(slave_seq)
                        else:
                            single_ctg = str(single_ctg) + 'NNNNN' + str(slave_seq)
                        SIN += 1
                        contam_size += len(slave_seq)
                        #for the log
                        gen_log.write(">" + str(slave) + '|' + str(slaveGene) + '|' + 'Single' +  "\n")
                        gen_log.write(str(slave_seq) + "\n")                       
                        #Increment
                        singelNumber += 1
                    if (singelNumberHGT < SingleHGT_size):
                        #Single HGT in sixth
                        #For the genome: Nothing to do
                        #slave
                        #Print gene name in the list for HGTsim and gene seq in a separate file for HGTsim
                        slave_ID = '_' + str(GC_cooDD[slave][slaveGene]['ID'])
                        slave_contig = slaveGene.replace(slave_ID, "")
                        slave_ctg = GC_genoDD[slave][slave_contig]
                        slave_start = int(GC_cooDD[slave][slaveGene]['start']) - 1
                        slave_end = int(GC_cooDD[slave][slaveGene]['end'])
                        slave_seq = slave_ctg[slave_start:slave_end]
                        slave_strand = GC_cooDD[slave][slaveGene]['strand']
                        if (slave_strand == '-1'):
                            seq = Seq(slave_seq)
                            reverse = seq.reverse_complement()
                            slave_seq = reverse
                        #For the genome
                        HGTsim_distribution.write(',' + str(slaveGene))
                        HGTsim_genes.write('>' + str(slaveGene) + "\n")
                        HGTsim_genes.write(str(slave_seq) + "\n")                       
                        SIN += 1
                        HGT_size += len(slave_seq)
                        #for the log
                        gen_log.write(">" + str(slave) + '|' + str(slaveGene) + '|' + 'SingleHGT-beforemutation' +  "\n")
                        gen_log.write(str(slave_seq) + "\n") 
                        #For HGT coor 
                        hgtcoo_log.write('chimeric-' + str(genumber) + "\t" + str(master) + "\t" +  str(slaveGene) 
                        + "\t" + 'Single' + "\t" + 'No-master-ctg' + "\t" + 'No-master-coo' + "\t" + 'No-master-coo' + "\n")                        
                        #Increment
                        singelNumberHGT += 1
                #Close the HGTsim distribution file
                HGTsim_distribution.write("\n")

                #Loop in genes of master, including Reeplacement, 
                #duplication and single are on sperate contig
                #Can be added sparately or at the end of the last master contig
                contig_of = GC_genoDD[master] 
                if (merge == 'no'):
                    #First print master contig, incld replacement
                    for contig in contig_of:
                        gen_out.write(">" + str(master) + '|' + str(contig) + "\n")
                        masterseq = contig_of[contig]
                        gen_out.write(str(masterseq) + "\n")
                        size = len(masterseq)
                        chim_size += size
                    if (Duplication_size > 0):
                        gen_out.write(">" + str(master) + '|' + 'DUPLICATE-CTG' + "\n")
                        gen_out.write(str(duplication_ctg) + "\n")
                        size = len(duplication_ctg)
                        chim_size += size
                    if (Single_size > 0):
                        gen_out.write(">" + str(master) + '|' + 'SINGLE-CTG' + "\n")
                        gen_out.write(str(single_ctg) + "\n")
                        size = len(duplication_ctg)
                        chim_size += size

                else:
                    #First print master contig, incld replacement
                    contigof_size = len(contig_of)
                    contignumber= 0
                    for contig in contig_of:
                        contignumber += 1
                        if (contignumber == contigof_size):
                            #last contig, merge duplicate and single to this contig
                            gen_out.write(">" + str(master) + '|' + str(contig) + "\n")
                            masterseq = contig_of[contig]
                            gen_out.write(str(masterseq))  
                            size = len(masterseq)
                            chim_size += size
                            if (Duplication_size > 0):
                                if (Single_size > 0):
                                    gen_out.write(str(duplication_ctg))
                                    size = len(duplication_ctg)
                                    chim_size += size
                                elif (Single_size == 0):
                                    gen_out.write(str(duplication_ctg) + "\n")
                                    size = len(duplication_ctg)
                                    chim_size += size
                            if (Single_size > 0):
                                gen_out.write(str(single_ctg) + "\n")
                                size = len(single_ctg)
                                chim_size += size
                        else:
                            #Not last contig, continue in printing contig
                            gen_out.write(">" + str(master) + '|' + str(contig) + "\n")
                            masterseq = contig_of[contig]
                            gen_out.write(str(masterseq) + "\n")
                            size = len(masterseq)
                            chim_size += size
                
                #Compute effective contams size
                totchim_size = chim_size + HGT_size
                totother_size = contam_size + HGT_size
                finalother = (float(totother_size) / float(totchim_size)) * 100
                finalhgt = (float(HGT_size) / float(totchim_size)) * 100
                finalcontam = (float(contam_size) / float(totchim_size)) * 100
                chim_log.write(str(master) + "\t" + str(taxoM) + "\t" + str(slave) + "\t" + str(taxoS)
                + "\t" + 'creation of chimeric genome: ID=' + str(genumber)
                + "\t" + 'Duplication=' + str(redundancy) + "\t" + 'Replacement=' + str(replacement) 
                + "\t" + 'Single=' + str(single) + "\t"
                + "\t" + 'DuplicationHGT=' + str(redundancyhgt) + "\t" + 'ReplacementHGT=' + str(replacementhgt)
                + "\t" + 'SingleHGT=' + str(singlehgt) 
                + "\t" + 'Chimeric_level=' + str(finalother) 
                + "\t" + 'Contamination_level=' + str(finalcontam) 
                + "\t" + 'HGT_level=' + str(finalhgt) + "\n")

    if (mode == 'HGT'):

        #Open mutant HGT genes file and store sequence
        mutant_file = open('input_sequence_mutant_nc.fasta')

        sequence = None
        ID = 'NA'
        mutant_of = {}

        for line in mutant_file:
            record = line.replace("\n", "")
            if '>' in line:
                if (ID != 'NA'):
                    #load seq
                    mutant_of[ID] = sequence
                    #reset
                    sequence = None
                    ID = 'NA'
                split_list = record.split(" ")
                ID = split_list[0]
                ID = ID.replace(">", "")
            else:
                #cat sequences
                if (sequence == None):
                    sequence = str(record)
                else:
                    sequence += str(record) 
        #Final seq load   
        mutant_of[ID] = sequence     


        #Open HGT coordinate file for replacement
        coo_file = open('HGT.coo')

        cooOthers_ofDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        cooReplace_ofDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        startReplace_ofDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        cor_ofDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        othertag = 0
        genomes_of = {}
        master_of = {}
        for line in coo_file:
            record = line.replace("\n", "")
            split_list = record.split("\t")
            genome = split_list[0]
            master = split_list[1]
            ID = split_list[2]
            category = split_list[3]
            contig = split_list[4]
            start = split_list[5]
            end = split_list[6]
            if (master != 'No-master'):
                master_of[genome] = master
            if (category == 'Replace'):
                cooReplace_ofDD[master][ID]['genome'] = genome
                cooReplace_ofDD[master][ID]['contig'] = contig
                cooReplace_ofDD[master][ID][contig]['start'] = start
                startReplace_ofDD[master][contig][ID] = start
                cor_ofDD[master][contig] = 0
                cooReplace_ofDD[master][ID][contig]['end'] = end
                cooReplace_ofDD[master][ID]['category'] = category
            else:
                cooOthers_ofDD[master][ID]['genome'] = genome
                cooOthers_ofDD[master][ID]['contig'] = contig
                cooOthers_ofDD[master][ID][contig]['start'] = start
                cooOthers_ofDD[master][ID][contig]['end'] = end
                cooOthers_ofDD[master][ID]['category'] = category
                othertag = 1
            genomes_of[genome] = 1

        #Open chimeric genomes file and strore sequences
        GC_genoDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        for genome in genomes_of:
            cat = str(genome) + '.fasta'
            genome_file = open(cat)
            defline = 'NA'
            sequence = None
            for line in genome_file:
                record = line.replace("\n", "")
                if '>' in line:
                    #fasta reading mode
                    if (defline != 'NA'):
                        #load the genes and reset
                        GC_genoDD[genome][defline] = sequence
                        defline = 'NA'
                        sequence = None
                    split_list = record.split(" ")
                    defline = split_list[0]
                    defline = defline.replace(">", "")
                    split_list = defline.split("|")
                    defline = split_list[1]
                else:
                    #cat sequences
                    if (sequence == None):
                        sequence = str(record)
                    else:
                        sequence += str(record)
            #final sequence load
            GC_genoDD[genome][defline] = sequence  

        #Proceed to replacement, first
        for master in startReplace_ofDD:
            start_ofDD = startReplace_ofDD[master]
            #Pass by contis
            for master_contig in start_ofDD:
                start_of = start_ofDD[master_contig]
                #Pass into sorted IDs
                for ID in sorted(start_of, key=start_of.get, reverse=False):
                    #Get mutant sequence
                    mutant_seq = mutant_of[ID]
                    mutantsize = len(mutant_seq)
                    #Get genome involved
                    genome = cooReplace_ofDD[master][ID]['genome']
                    #Get sequence of 
                    sequence = GC_genoDD[genome][master_contig]
                    #Get coordinate of original sequence
                    master_cor = cor_ofDD[master][master_contig]
                    Ostart = int(cooReplace_ofDD[master][ID][master_contig]['start']) + int(master_cor)
                    Oend = int(cooReplace_ofDD[master][ID][master_contig]['end']) + int(master_cor)
                    #Get Original sequence
                    Oseq = sequence[Ostart:Oend]
                    Osize = len(Oseq)
                    #Compuet size difference a modify cor factor
                    mod_size = mutantsize - Osize
                    mod_cor = int(master_cor) + int(mod_size)
                    cor_ofDD[master][master_contig] = mod_cor
                    #Insert mutant sequence
                    chim_seq = sequence[:Ostart] + str(mutant_seq) + sequence[Oend:]
                    #Reload chim seq in hash
                    GC_genoDD[genome][master_contig] = chim_seq

        #Proceed to duplication
        SUPP_genoDD = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        #First check if any duplication or single activated
        if (othertag == 1): 
            for master in cooOthers_ofDD:
                supplemental_ctg = None
                ID_ofDD = cooOthers_ofDD[master]
                for ID in ID_ofDD:
                    #Get mutant sequence
                    mutant_seq = mutant_of[ID]
                    #If not replacement, add it into supp contig
                    if (supplemental_ctg == None):
                        supplemental_ctg = 'NNNNN' + str(mutant_seq)
                    else:
                        supplemental_ctg = str(supplemental_ctg) + 'NNNNN' + str(mutant_seq)
                    SUPP_genoDD[genome]['SUPP'] = supplemental_ctg 

        #Print new chimeric genome
        for genome in genomes_of:
            #Open outfile
            cat = str(genome) + '-HGT.ali'
            gen_out = open(cat, "w")

            #Get info / seq
            master = master_of[genome]
            contig_of = GC_genoDD[genome]
            supp_of = {}
            if (genome in SUPP_genoDD):
                supp_of = SUPP_genoDD[genome]
            #Pass into contig of master / supp contig and print
            if (merge == 'no'):
                #First print master contig, incld replacement
                for contig in contig_of:
                    gen_out.write(">" + str(master) + '|' + str(contig) + "\n")
                    masterseq = contig_of[contig]
                    gen_out.write(str(masterseq) + "\n")
                if ('SUPP' in supp_of):
                    supp_ctg = supp_of['SUPP']
                    gen_out.write(">" + str(master) + '|' + 'SUPP-CTG' + "\n")
                    gen_out.write(str(supp_ctg) + "\n")
            else:
                #First print master contig, incld replacement
                contigof_size = len(contig_of)
                contignumber= 0
                for contig in contig_of:
                    contignumber += 1
                    if (contignumber == contigof_size):
                        #last contig, merge duplicate and single to this contig
                        gen_out.write(">" + str(master) + '|' + str(contig) + "\n")
                        masterseq = contig_of[contig]
                        gen_out.write(str(masterseq))  
                        #Print supp contig in merge
                        if ('SUPP' in supp_of):
                            supp_ctg = supp_of['SUPP']
                            gen_out.write(str(supp_ctg) + "\n")
                    else:
                        #Not last contig, continue in printing contig
                        gen_out.write(">" + str(master) + '|' + str(contig) + "\n")
                        masterseq = contig_of[contig]
                        gen_out.write(str(masterseq) + "\n")
           


if __name__ == '__main__':
    main()

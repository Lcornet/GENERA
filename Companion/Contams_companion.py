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
from collections import defaultdict

@click.command()
###ARGUMENT#####
#main file - depending of the mode can be checkm result or rnammer output
@click.argument('main_file', type=click.Path(exists=True,readable=True))
#taxid
@click.option('--mode', default='no', help='mode')
@click.option('--ckcompleteness', default='95', help='Checkm Comp')
@click.option('--ckcontamination', default='5', help='Checkm Comp')
@click.option('--gunccss', default='0.1', help='GUNC css')
@click.option('--guncrrs', default='0.5', help='GUNC rrs')
@click.option('--physetercontamination', default='100', help='Physeter conta')
@click.option('--krakencontamination', default='100', help='Kraken conta')
@click.option('--bucompleteness', default='0', help='Physeter conta')
@click.option('--budups', default='100', help='Kraken conta')
@click.option('--ck2completeness', default='95', help='Checkm2 compl')
@click.option('--ck2contamination', default='5', help='Checkm2 conta')
@click.option('--numcontigs', default='100', help='Number of contigs')

def main(main_file, mode, ckcompleteness, ckcontamination, gunccss, guncrrs, physetercontamination, krakencontamination, bucompleteness, budups, ck2completeness, ck2contamination, numcontigs):
		
	#read fasta file
    if (mode == 'table'):
        file_out = open('GENERA-contams.table', "w")
        list_out = open('positive-list.txt', "w")
        contam_of = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        compl_of = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))
        quast_of = defaultdict( lambda: defaultdict(lambda: defaultdict( dict )))

        #Open checkm results
        CKfile = open('results.Checkm')
        for line in CKfile:
            record = line.replace("\n", "")
            if ('#' in record):
                continue
            else:
                split_list = record.split(",")
                Genome = split_list[0]
                completeness = split_list[1]
                contam = split_list[2]
                Str = split_list[3]
                contam_of[Genome]['checkm'] = contam
                contam_of[Genome]['checkmStr'] = Str
                compl_of[Genome]['checkm'] = completeness
        

        #Open gunc results
        GUfile = open('GUNC.progenomes_2.1.maxCSS_level.tsv')
        for line in GUfile:
            record = line.replace("\n", "")
            if ('n_genes_called' in record):
                continue
            else:
                split_list = record.split("\t")
                Genome = split_list[0]
                CSS = split_list[7]
                contam= split_list[8]
                RRS = split_list[11]
                statu = split_list[12]
                contam_of[Genome]['gunc'] = contam
                contam_of[Genome]['gunc-RRS'] = RRS
                contam_of[Genome]['gunc-CSS'] = CSS 
                contam_of[Genome]['gunc-status'] = statu 

        #Open busco results
        BUfile = open('batch_summary.txt')
        for line in BUfile:
            record = line.replace("\n", "")
            if (('Input_file' in record) or ('logs' in record)):
                continue
            else:
                split_list = record.split("\t")
                Genome = split_list[0]
                Genome = Genome.replace(".fna", "")
                placement = split_list[1]
                compl = split_list[2]
                dups = split_list[4]
                contam_of[Genome]['busco'] = dups
                contam_of[Genome]['busco-place'] = placement
                compl_of[Genome]['busco'] = compl

        #Physeter results
        PYfile = open('Physeter.report')
        for line in PYfile:
            record = line.replace("\n", "")
            split_list = record.split("\t")
            Genome = split_list[0]
            Genome = Genome.replace("-split-kraken", "")
            placement = split_list[1]
            conta = split_list[3]
            contam_of[Genome]['physeter'] = conta
            contam_of[Genome]['physeter-place'] = placement

        #Kraken results
        KAfile = open('Kraken.report')
        for line in KAfile:
            record = line.replace("\n", "")
            split_list = record.split("\t")
            Genome = split_list[0]
            Genome = Genome.replace("-split", "")
            placement = split_list[1]
            conta = split_list[3]
            contam_of[Genome]['kraken'] = conta
            contam_of[Genome]['kraken-place'] = placement       
    

        #Quast results
        QUfile = open('quast.report')
        for line in QUfile:
            record = line.replace("\n", "")
            if ('Assembly' in record):
                continue
            else:
                split_list = record.split("\t")
                Genome = split_list[0]
                contigs = split_list[1]
                length = split_list[7]
                GC = split_list[16]
                N50 = split_list[17]
                quast_of[Genome]['contigs'] = contigs
                quast_of[Genome]['length'] = length
                quast_of[Genome]['GC'] = GC
                quast_of[Genome]['N50'] = N50
        
        #Checkm2
        CK2file = open('quality_report.tsv')
        for line in CK2file:
            record = line.replace("\n", "")
            if ('Name' in record):
                continue
            else:
                split_list = record.split("\t")
                Genome = split_list[0]
                completeness = split_list[1]
                conta = split_list[2]
                contam_of[Genome]['CK2'] = conta
                compl_of[Genome]['CK2'] = completeness  

        
        #print(quast_of)

        #Write the table
        file_out.write('Genome' + "\t" + 'Checkm_completeness' + "\t" + 'Checkm_contamination' + "\t" + 'Checkm_str' + "\t" 
        + 'Busco_placemenent' + "\t" + 'Busco_completeness' + "\t" + 'Busco_duplicate' 
        + "\t" + 'GUNC_CSS' + "\t" + 'GUNC_RRS' + "\t" + 'GUNC_status'
        + "\t" + 'GUNC_conta' + "\t" + 'Physeter_placement' + "\t" + 'Physeter_contamination' 
        + "\t" + 'Kraken_placement' + "\t" + 'Kraken_contamination'
        + "\t" + 'Checkm2_completeness' + "\t" + 'Checkm2_contamination'
        + "\t" + "quast_#contigs" + "\t" + 'quast_tot_length' + "\t" + 'quast_GC' + "\t" + 'quast_N50' + "\n")     

        for genome in contam_of:
            #print(genome)
            tool_of = contam_of[genome]
            Checkm_comp = compl_of[genome]['checkm']
            Busco_comp = compl_of[genome]['busco']
            Checkm_conta = tool_of['checkm']
            Checkm_str = tool_of['checkmStr']
            Busco_conta = tool_of['busco']
            Busco_place = tool_of['busco-place']
            GUNC_conta = tool_of['gunc']
            GUNC_RRS = tool_of['gunc-RRS']
            GUNC_CSS = tool_of['gunc-CSS']
            GUNC_status = tool_of['gunc-status']
            Physeter_conta = tool_of['physeter']
            Physeter_place = tool_of['physeter-place']
            Kraken_conta = tool_of['kraken']
            Kraken_place = tool_of['kraken-place']
            Checkm2_conta = tool_of['CK2']
            Checkm2_comp = compl_of[genome]['CK2']
            quast_contigs = quast_of[genome]['contigs']
            quast_length = quast_of[genome]['length']
            quast_GC = quast_of[genome]['GC']
            quast_N50 = quast_of[genome]['N50']

            file_out.write(str(genome) + "\t" + str(Checkm_comp) + "\t" + str(Checkm_conta) + "\t" + str(Checkm_str) 
            + "\t" + str(Busco_place) + "\t" + str(Busco_comp) + "\t" + str(Busco_conta) 
            + "\t" + str(GUNC_CSS) + "\t" + str(GUNC_RRS) + "\t" + str(GUNC_status) 
            + "\t" + str(GUNC_conta) + "\t" + str(Physeter_place) + "\t" + str(Physeter_conta) 
            + "\t" + str(Kraken_place) + "\t" + str(Kraken_conta)
            + "\t" + str(Checkm2_comp) + "\t" + str(Checkm2_conta)
            + "\t" + str(quast_contigs) + "\t" + str(quast_length) + "\t" + str(quast_GC) + "\t" + str(quast_N50) + "\n")

            #print in positive list
            #print(contam_of)
            if ((float(Checkm_comp) > float(ckcompleteness)) and (float(Checkm_conta) < float(ckcontamination))
            and (float(Busco_comp) > float(bucompleteness)) and (float(Busco_conta) < float(budups))
            and (float(GUNC_CSS) < float(gunccss)) and (float(GUNC_RRS) > float(guncrrs))
            and (int(quast_contigs) < float(numcontigs)) 
            and (float(Physeter_conta) < float(physetercontamination)) 
            and (float(Kraken_conta) < float(krakencontamination))
            and (float(Checkm2_conta) < float(ck2contamination))
            and (float(Checkm2_comp) > float(ck2completeness))):
                list_out.write(str(genome) + "\n")

if __name__ == '__main__':
    main()
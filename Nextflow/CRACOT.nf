#!/usr/bin/env nextflow

/*
========================================================================================
                        GENERA contams
========================================================================================
GIT url : https://github.com/Lcornet/GENERA
----------------------------------------------------------------------------------------
*/


/*
HELP SECTION
*/

//Define your help message, this section is supposed
//to state how to use the workflow
//Can be as long as needed, no tabulation, only white-space here
def helpMessage() {
	log.info """

    Description:

    Version: 1.0.0 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow run CRACOT.nf --genomes=genomes --lineage=Genomes.taxomonomy --list=positive-list.txt 
    --cpu=60 --num=100 --taxorank=phylum --redundant=2 --replacement=2 --single=4 
    --hgtrate=none --hgtrandom=no --redundanthgt=0 --replacementhgt=0 --singlehgt=0 
    
    Mandatory arguments:
    --genomes                Specify directory with genomes
    --lineage                Specify lineage file path
    --list                   Specify positive list file path

    Optional arguments:
    --taxorank               Specify taxonomic rank, phylum or order or family or genus or species, default=phylum
    --num                    Specify number of chimeric genomes, default=100
    --maskslave              Used only slave genome one time, default = no
    --redundant              Specify number of redundant events, default=5
    --replacement            Specify number of replacement events, default=5
    --single                 Specify number of single events, default=5
    --hgtrate                Specify the mutation rate (1 to 99%) of HGT events, default = none
    --redundanthgt           Specify number of hgt redundant events, default=0
    --replacementhgt         Specify number of hgt replacement events, default=0
    --singlehgt              Specify number of hgt single events, default=0   
    --merge                  Merge redundant and single to the last contig of the chimeric genome, yes or no, default = yes
    --hgtrandom              Activate the random insertion of HGT events, yes or no, default = no, Recommended with SINGLE hgt events only.
    --cpu                    number of cpus to use, default = 1

    """.stripIndent()
    
}
 
// Show help message and exit the workflow (--help option)
params.help = null
if (params.help){
    helpMessage()
    exit 0
}

/*
INPUT AND OPTIONS SETTING
*/

//Path to OGs : Mandatory
params.genomes = null
if (params.genomes == null) {
	exit 1, "Specify directory with genomes"
}

//Path to OGs : Mandatory
params.lineage = null
if (params.lineage == null) {
	exit 1, "Specify lineage file path"
}

//Path to OGs : Mandatory
params.list = null
if (params.list == null) {
	exit 1, "Specify positive list file path"
}

//taxolevel
params.taxorank = 'phylum'

//num
params.num = '100'

//maskslave 
params.maskslave = 'no'

//duplication
params.redundant = '5'

//replacement
params.replacement = '5'

//single
params.single = '5'

//merge
params.merge = 'yes'

//activate HGT
params.hgtrate = 'none'

//duplication HGT
params.redundanthgt = '0'

//replacement HGT
params.replacementhgt = '0'

//single HGT
params.singlehgt = '0'

//HGT random insertion
params.hgtrandom = 'no'

//cpu
params.cpu = '1'

//outdir
params.outdir='GENERA-chimeric'

//Path to companion
params.companion = '/opt/companion/chimeric_companion.py'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
genomes_ch1 = Channel.fromPath(params.genomes)
genomes_ch2 = Channel.fromPath(params.genomes)
lineage_ch = Channel.fromPath(params.lineage)
list_ch = Channel.fromPath(params.list)

//list
process lineage {
	//informations

	//input output
    input:
    file 'list' from list_ch
    file 'lineage' from lineage_ch

    output:
    file 'lineageSub' into sublineage_ch
    file 'list' into list_ch2


    //script
    script:
    """
    for f in `cat list`; do grep \$f lineage; done > lineageSub
    """
}

//makecorr
process makecorr {
	//informations

	//input output
    input:
    file 'lineageSub' from sublineage_ch
    val num from params.num
    val companion from params.companion
    val level from params.taxorank

    output:
    file 'chimeric.corr' into corr_ch
    file 'chimeric.list' into toof_ch1
    file 'chimeric.list' into toof_ch2
    file 'chimeric.list' into toof_ch3
    file 'chimeric.idl' into idl_ch

    //script
    script:
    if (params.maskslave == 'yes'){ 
        """
        $companion  lineageSub --mode=list --level=$level --number=$num --mask_slave=yes
        """
    }
    else if (params.maskslave == 'no'){
        """
        $companion  lineageSub --mode=list --level=$level --number=$num --mask_slave=no
        """
    }
}

//plasmiddel
process plasmiddel {
	//informations

	//input output
    input:
    file '*' from genomes_ch1
    file 'chimeric.list' from toof_ch1

    output:
    file 'SUB' into subgenome_ch1
    file 'SUB' into subgenome_ch2
    

    //script
    script:
    """
    mkdir PLA
    mkdir SUB
    for f in `cat chimeric.list`; do cp genomes/\$f.fna .; done
    #for f in `cat chimeric.list`; do head -n5 \$f.fna > \$f.head; done
    #for f in `cat chimeric.list`; do mkdir plas; cp \$f.fna plas; plasmidpicker pick -i \$f.fna -d plas; done
    for f in `cat chimeric.list`; do mkdir plas; cp \$f.fna plas; plasmidpicker pick -i \$f.fna -d plas; \
    grep ">" plas/plasmids*.fna | cut -f2 -d"|" | cut -f1 -d" " > plasmid.id; \
    fasta2ali.pl \$f.fna; grep ">" \$f.ali > all-id.list; \
    for d in `cat plasmid.id`; do grep -v \$d all-id.list > temp; mv temp all-id.list; done; \
    for e in `cat all-id.list`; do grep -A1 \$e \$f.ali; done > \$f-sub.ali; ali2fasta.pl \$f-sub.ali; \
    mv \$f-sub.fasta SUB/\$f.fna; rm -rf plas plasmid.id all-id.list; done
    """

}


//prodigal
process prodigal {
	//informations

	//input output
    input:
    //file '*' from genomes_ch1
    file 'SUB' from subgenome_ch1
    file 'chimeric.list' from toof_ch2

    output:
    file 'FAA/' into faa_ch
    file 'GENES/' into genes_ch

    //script
    script:
    """
    mkdir GEN
    #cp genomes/*.fna GEN/
    cp SUB/*.fna GEN/
    for f in `cat chimeric.list`; do prodigal -i GEN/\$f.fna -o \$f.genes -d \$f.genes.fna -a \$f.faa; done
    mkdir GENES
    mkdir FAA
    mv *.genes.fna GENES/
    mv *.faa FAA/
    """

}

//of
process orthofinder {
	//information

	//input output
    input:
    file '*' from faa_ch
    val cpu from params.cpu

    output:
    file 'OG/' into og_ch

    //script
    script:
    """
    mkdir OF-indir
    cp FAA/*.faa OF-indir/
    orthofinder.py -t $cpu -a $cpu -f OF-indir
    mkdir OG
    cp OF-indir/OrthoFinder/Results_*/Orthogroup_Sequences/*.fa OG/
    """
}


//makechim
process makechim {
	//informations

	//input output
    input:
    file '*' from genomes_ch2
    file '*' from og_ch
    file '*' from genes_ch
    file 'chimeric.corr' from corr_ch
    file 'chimeric.list' from toof_ch3
    val duplication from params.redundant
    val replacement from params.replacement
    val single from params.single
    val merge from params.merge
    val duplicationhgt from params.redundanthgt
    val replacementhgt from params.replacementhgt
    val singlehgt from params.singlehgt
    val companion from params.companion

    output:
    file 'chimeric.log' into chimlog_ch
    file 'CHIMERIC-genomes/' into chimgen_ch
    file 'CHIMERIC-sequences/' into chimgenlog_ch
    file 'CHIMERIC-hgtsim/' into chimhgt_ch
    file 'all-HGT.coo' into hgtcoo_ch

    //script
    script:
    """
    for f in `cat chimeric.list`; do cp genomes/\$f*.fna .; done
    cp chimeric.corr chimeric.corr.bak
    cp GENES/*.fna .
    $companion chimeric.list --mode=chimeric --corr_file=chimeric.corr --redundancy=$duplication --replacement=$replacement \
    --single=$single --merge=$merge --redundancyhgt=$duplicationhgt --replacementhgt=$replacementhgt --singlehgt=$singlehgt
    mkdir CHIMERIC-genomes/
    mkdir CHIMERIC-sequences/
    mkdir CHIMERIC-hgtsim/
    echo 'RUN companion in chimeric mode' > chimericINFO-log.ali
    echo 'RUN companion in chimeric mode' > chimericINFO.ali
    mv chimeric*-log.ali CHIMERIC-sequences/
    mv chimeric-*-genes.ali CHIMERIC-hgtsim/
    mv chimeric*.ali CHIMERIC-genomes/
    mv chimeric-*-distribution.txt CHIMERIC-hgtsim/
    #HGT coo file
    cat *HGT.coo > all-HGT.coo
    """  
}

//hgtsim
process hgtsim {
	//informations

	//input output
    input:
    file 'chimeric.log' from chimlog_ch
    file 'CHIMERIC-genomes/' from chimgen_ch
    file 'CHIMERIC-sequences/' from chimgenlog_ch 
    file 'CHIMERIC-hgtsim/' from chimhgt_ch
    file 'chimeric.idl' from idl_ch 
    val rat from params.hgtrate

    output:
    file 'CHIMERIC/' into chimF_ch
    file 'CHIM-log/' into chimseqF_ch
    file 'CHIM-hgt/' into chimhgtF_ch
    file 'chimeric-genomes.list' into chimgenF_ch
    file 'chimeric-use.idl' into idlF_ch
    file 'input_sequence_mutant_nc.fasta' into mutantseq_ch
    file 'HGTSIM-GENOMES' into hgtsimF_ch

    //script
    script:
    if (params.hgtrate == 'none') {
        """
        mkdir CHIMERIC
        mkdir CHIM-log
        cp -r CHIMERIC-genomes/ CHIMERIC-genomes-R/
        cd CHIMERIC-genomes-R/CHIMERIC-genomes/
        rm -f *INFO*
        ali2fasta.pl *.ali
        cd ../../
        mv CHIMERIC-genomes-R/CHIMERIC-genomes/*.fasta CHIMERIC/
        cp -r CHIMERIC-sequences/ CHIMERIC-sequences-R/
        cd CHIMERIC-sequences-R/CHIMERIC-sequences/
        rm -f *INFO*
        ali2fasta.pl *.ali 
        cd ../../
        mv CHIMERIC-sequences-R/CHIMERIC-sequences/*.fasta CHIM-log/
        mv chimeric.log chimeric-genomes.list
        mv chimeric.idl chimeric-use.idl
        #No HGT, just store files
        mkdir CHIM-hgt
        mv CHIMERIC-hgtsim/* CHIM-hgt/
        mkdir HGTSIM-GENOMES
        echo 'NO HGT' > HGTSIM-GENOMES/info.txt
        echo 'NO HGT' > input_sequence_mutant_nc.fasta
        """       
    }
    else {
        """
        #genes part
        cat CHIMERIC-hgtsim/CHIMERIC-hgtsim/*genes.ali > genes-all.ali
        mv genes-all.ali genes.ali
        ali2fasta.pl genes.ali  

        #Distribution file part
        cat CHIMERIC-hgtsim/CHIMERIC-hgtsim/chimeric-*-distribution.txt > distribution.txt

        #Genomes dir part
        mkdir input_genomes
        cp -r CHIMERIC-genomes/ CHIMERIC-genomes-R/
        cd CHIMERIC-genomes-R/CHIMERIC-genomes/
        rm -f *INFO*
        ali2fasta.pl *.ali
        cd ../../
        cp CHIMERIC-genomes-R/CHIMERIC-genomes/*.fasta input_genomes/
        mkdir CHIMERIC
        mv CHIMERIC-genomes-R/CHIMERIC-genomes/*.fasta CHIMERIC/
   

        #Run HGT sim
        mkdir HGTSIM-GENOMES
        echo HgtSIM_outputs_$rat\_1-0-1-1 > dir
        mkdir \$(<dir)
        #mkdir HgtSIM_outputs_10_1-0-1-1/
        HgtSIM -t genes.fasta -d distribution.txt -f input_genomes -r 1-0-1-1 -x fasta -i $rat 
        mv \$(<dir)/Genomes_with_transfers/* HGTSIM-GENOMES/

        #log sequence part
        mkdir CHIM-log
        cp -r CHIMERIC-sequences/ CHIMERIC-sequences-R/
        cd CHIMERIC-sequences-R/CHIMERIC-sequences/
        rm -f *INFO*
        ali2fasta.pl *.ali 
        cd ../../
        mv CHIMERIC-sequences-R/CHIMERIC-sequences/*.fasta CHIM-log/

        #log part
        mv chimeric.log chimeric-genomes.list
        mv chimeric.idl chimeric-use.idl
        mkdir CHIM-hgt
        mv \$(<dir)/*.txt CHIM-hgt/
        #mutant seq
        cp \$(<dir)/input_sequence_mutant_nc.fasta .
        """  
    }
}

//HGT introduce
process hgtintroduce {
	//information

	//input output
    input:
    file 'CHIMERIC/' from chimF_ch
    file 'CHIM-log/' from chimseqF_ch
    file 'CHIM-hgt/' from chimhgtF_ch
    file 'chimeric-genome.list' from chimgenF_ch
    file 'chimeric-use.idl' from idlF_ch 
    file 'all-HGT.coo' from hgtcoo_ch
    file 'input_sequence_mutant_nc.fasta' from mutantseq_ch
    file 'HGTSIM-GENOMES' from hgtsimF_ch
    val companion from params.companion
    val merge from params.merge

    output:
    file 'FINAL-genomes/' into chimI_ch
    file 'FINAL-sequences/' into chimseqI_ch
    file 'FINAL-info/' into chimhgtI_ch
    file 'chimeric-genomes-FINAL.list' into chimgenI_ch
    file 'chimeric-FINAL.idl' into idlI_ch 

    //script
    script:
    if (params.hgtrate == 'none') {
        """
        #NO HGT to treat, just copy genomes
        mkdir FINAL-genomes/
        cp CHIMERIC/CHIMERIC/* FINAL-genomes/
        mkdir FINAL-sequences
        cp CHIM-log/CHIM-log/* FINAL-sequences/
        mkdir FINAL-info/
        cp CHIM-hgt/CHIM-hgt/CHIMERIC-hgtsim/* FINAL-info/
        cp -r chimeric-genome.list chimeric-genomes-FINAL.list
        cp -r chimeric-use.idl chimeric-FINAL.idl
        """
    }
    else if (params.hgtrandom == 'yes') {
        """
        #Activation of the radnom insertion from HGTsim
        #Copy files from HGTSIM genomes dir
        mkdir FINAL-genomes/
        mv HGTSIM-GENOMES/*.fasta FINAL-genomes/
        mkdir FINAL-sequences
        cp CHIM-log/CHIM-log/* FINAL-sequences/
        mkdir FINAL-info/
        cp CHIM-hgt/CHIM-hgt/* FINAL-info/
        cp -r chimeric-genome.list chimeric-genomes-FINAL.list
        cp -r chimeric-use.idl chimeric-FINAL.idl
        """
    }
    else {
        """
        cp all-HGT.coo HGT.coo
        cp CHIMERIC/CHIMERIC/* .
        $companion HGT.coo --mode=HGT --merge=$merge
        rm -f *.fasta
        ali2fasta.pl *HGT.ali
        find *-HGT.fasta > list
        sed -i -e 's/-HGT.fasta//g' list
        for f in `cat list`; do mv \$f-HGT.fasta \$f.fasta; done
        #Prepare copy
        mkdir FINAL-genomes/
        mv *.fasta FINAL-genomes/
        mkdir FINAL-sequences
        cp CHIM-log/CHIM-log/* FINAL-sequences/
        mkdir FINAL-info/
        cp CHIM-hgt/CHIM-hgt/* FINAL-info/
        cp -r chimeric-genome.list chimeric-genomes-FINAL.list
        cp -r chimeric-use.idl chimeric-FINAL.idl
        """
    }
}

//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'FINAL-genomes/' from chimI_ch
    file 'FINAL-sequences/' from chimseqI_ch
    file 'FINAL-info/' from chimhgtI_ch
    file 'chimeric-genomes-FINAL.list' from chimgenI_ch
    file 'chimeric-FINAL.idl' from idlI_ch 
    file 'SUB' from subgenome_ch2

    output:
    file 'CHIMERIC-genomes/' into chimFinal_ch
    file 'CHIM-sequences/' into chimseqFinal_ch
    file 'CHIM-hgt-info/' into chimhgtinfoFinal_ch
    file 'GENOMES-USED_out-plasmid/' into genomesusedFinal_ch
    file 'chimeric-genomes.list' into chimgenFinal_ch
    file 'chimeric.idl' into idlFinal_ch

    //script
    script:
    """
    mkdir CHIMERIC-genomes/
    cp FINAL-genomes/FINAL-genomes/* CHIMERIC-genomes/
    mkdir CHIM-sequences/
    cp FINAL-sequences/FINAL-sequences/* CHIM-sequences/
    mkdir CHIM-hgt-info/
    cp FINAL-info/FINAL-info/* CHIM-hgt-info/
    cp -r chimeric-genomes-FINAL.list chimeric-genomes.list
    cp -r chimeric-FINAL.idl chimeric.idl
    mkdir GENOMES-USED_out-plasmid
    cp SUB/* GENOMES-USED_out-plasmid/
    """
}
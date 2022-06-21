#!/usr/bin/env nextflow

/*
========================================================================================
                         Phylogeny GENERA
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

    Version: 2.0.0 
    
    Citation:
    Please cite : 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow Phylogeny.nf --OG=OGs --IDM=file.idm --jackk=yes 
    
    Mandatory arguments:
    --OG                    Path to OG directory in fasta format (.faa for prot and .fna files for DNA)
    --IDM                   Path to IDM file

    Optional arguments:
    --mode                  specify prot or DNA, default = prot
    --degap                 activate degap of sequences (delete all gap befaore new alignment), yes or no, default = no
    --align                 activate the alignment for protein, yes or no, default = no
    --cpu                   number of cpus to use, default = 1

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
params.OG = null
if (params.OG == null) {
	exit 1, "Path to OG directory in fasta format."
}

//Path to IDM
params.IDM = null
if (params.IDM == null) {
	exit 1, "Path to IDM file containing final names."
}

//mode
params.mode = 'prot'

//degap
params.degap = 'no'

//align
params.align = 'no'

//outdir
params.outdir='GENERA-phylogeny-single'

//cpu
params.cpu = '1'

//Path to companion
params.companion = '/opt/companion/Phylogeny_companion.py'

//version
params.version = '2.0.0'

/*
CORE PROGRAM
*/

//Load input files
og_ch = Channel.fromPath(params.OG)
idm_ch = Channel.fromPath(params.IDM)


//BMC names format
process format {
	//informations

	//input output
    input:
    file '*.f*' from og_ch

    output:
    file 'FASTA/*abbr*' into abbr_ch1
    file "GENERA-Phylogeny.log" into log_ch0

    //script
    script:
    println "GENERA info: format names with BMC"
    if (params.mode == 'prot') {
        println "GENERA info: Phylogeny protein mode"
        """
        mkdir FASTA
        cp .f/* FASTA/
        cd FASTA
        find *.* | cut -f1 -d'.' > list
        for f in `cat list`; do sed -i -e 's/@/ /g' \$f*.faa; inst-abbr-ids.pl \$f*.faa --id-regex=:DEF; sed -i -e 's/|//g' \$f*abbr*.faa; done
        cd ../
        echo "GENERA info: Phylogeny protein mode" >> GENERA-Phylogeny.log
        echo "GENERA info: format names with BMC" >> GENERA-Phylogeny.log
        """
    }
    else if (params.mode == 'DNA') {
        println "GENERA info: Phylogeny DNA mode"
        """
        mkdir FASTA
        cp .f/* FASTA/
        cd FASTA
        find *.* | cut -f1 -d'.' > list
        for f in `cat list`; do sed -i -e 's/@/ /g' \$f*.fna; inst-abbr-ids.pl \$f*.fna --id-regex=:DEF; sed -i -e 's/|//g' \$f*abbr*.fna; done
        cd ../
        echo "GENERA info: Phylogeny DNA mode" >> GENERA-Phylogeny.log
        echo "GENERA info: format names with BMC" >> GENERA-Phylogeny.log
        """      
    }
}


//degap
process degap {
	//informations

	//input output
    input:
    file '*' from abbr_ch1
    file "GENERA-Phylogeny.log" from log_ch0

    output:
    file 'DEGAP' into degap_ch
    file "GENERA-Phylogeny.log" into log_ch1

    //script
    script:
    println "GENERA info: degap part"
    if (params.degap == 'yes') {
        if (params.mode == 'prot') {
            """
            find *.faa | cut -f1 -d"." > list
            for f in `cat list`; do fasta2ali.pl \$f.faa; rm -f \$f.faa; done
            for f in `cat list`; do ali2fasta.pl \$f.ali --degap; done
            for f in `cat list`; do fasta2ali.pl \$f.fasta; mv \$f.fasta \$f.faa; done
            mkdir DEGAP
            mv *.faa DEGAP/
            echo "GENERA info: run degap mode" >> GENERA-Phylogeny.log
            """
        }
        else if (params.mode == 'DNA') {
            """
            find *.fna | cut -f1 -d"." > list
            for f in `cat list`; do fasta2ali.pl \$f.fna; rm -f \$f.fna; done
            for f in `cat list`; do ali2fasta.pl \$f.ali --degap; done
            for f in `cat list`; do fasta2ali.pl \$f.fasta; mv \$f.fasta \$f.fna; done
            mkdir DEGAP
            mv *.fna DEGAP/
            echo "GENERA info: run degap mode" >> GENERA-Phylogeny.log
            """
        }
    }
    else if (params.degap == 'no') {
        """
        mkdir DEGAP
        mv *.f* DEGAP/
        echo "GENERA info: degap mode not activated" >> GENERA-Phylogeny.log
        """            
    }
}

//Alignment
process alignment {
	//informations

	//input output
    input:
    file 'DEGAP' from degap_ch
    file "GENERA-Orthology.log" from log_ch1

    output:
    file 'aligned/*' into ali_ch1
    file "GENERA-Phylogeny.log" into log_ch2

    //script
    script:
    println "GENERA info: run alignment"
    if (params.mode == 'prot') {
        if (params.align == 'yes') {
            """
            mv DEGAP/* .
            mkdir aligned
            mkdir OGs
            find *.faa | cut -f1 -d"." > list
            mv *.faa OGs/
            for f in `cat list`; do muscle3.8.31_i86linux64 -in OGs/\$f.faa -out aligned/\$f.faa; done
            echo "GENERA info: run prot alignment" >> GENERA-Phylogeny.log
            """
        }
        else if (params.align == 'no') {
            """
            mv DEGAP/* .
            mkdir aligned
            mv *.faa aligned/
            echo "GENERA info: No prot alignment" >> GENERA-Phylogeny.log
            """            
        }
    }
    else if (params.mode == 'DNA') {
        """
        mv DEGAP/* .
        mkdir aligned
        mv *.fna aligned/
        echo "GENERA info: No alignment for DNA files" >> GENERA-Phylogeny.log
        """
    }
}

//Unambiguous position selection
process unambiguousPosition {
	//informations

	//input output
    input:
    file '*' from ali_ch1
    file "GENERA-Phylogeny.log" from log_ch2

    output:
    file 'UNAMBIG/*.ali' into unambigous_ch1
    file 'UNAMBIG/*.ali' into unambigous_ch2
    file 'ALIGNU' into alignu_ch
    file "GENERA-Phylogeny.log" into log_ch3

    //script
    script:
    println "GENERA info: Unambiguous position conservation"
    if (params.mode == 'prot') {
        """
        mkdir BMGE
        mkdir UNAMBIG
        mkdir ALIGNU
        mv *.faa BMGE/
        cd BMGE
        fasta2ali.pl *.faa
        ali2phylip.pl *.ali --bmge-mask=medium --ali
        ali2fasta.pl *.ali
        cd ../
        mv BMGE/*a2p.ali UNAMBIG/
        cp BMGE/*a2p.fasta ALIGNU/
        echo "GENERA info: Unambiguous position conservation" >> GENERA-Phylogeny.log
        """
    }
    else if (params.mode == 'DNA') {
        """
        mkdir UNAMBIG
        mkdir ALIGNU
        mv *.fna UNAMBIG
        cd UNAMBIG/
        fasta2ali.pl *.fna
        cd ../
        cp *.fna ALIGNU/
        echo "GENERA info: No unambiguous position conservation for DNA files" >> GENERA-Phylogeny.log
        """
    }
}

//ML inference
process mlinference {
	//informations

	//input output
    input:
    file '*' from unambigous_ch1
    val cpu from params.cpu
    file "GENERA-Phylogeny.log" from log_ch3

    output:
    file 'TREE' into mltree_ch1
    file "GENERA-Phylogeny.log" into log_ch4

    //script
    script:
    if (params.mode == 'prot') {
        println "GENERA info: Prot ML inference"
        """
        mkdir ali
        mv *.ali ali/
        cd ali/
        find *.ali > list
        sed -i -e 's/.ali//g' list
        for f in `cat list`; do ali2phylip.pl \$f.ali --map-ids; done
        for f in `cat list`; do raxmlHPC-PTHREADS-AVX -T $cpu -s \$f.phy \
        -n \$f-RAXML-PROTGAMMALGF-100xRAPIDBP -m PROTGAMMALGF -N 100 -f a -x 1975021703574 -p 1975021703574; done
        #format
        for f in `cat list`; do mv \$f.idm RAxML_bipartitions.idm; \
        format-tree.pl RAxML_bipartitions.\$f-RAXML-PROTGAMMALGF-100xRAPIDBP --map-ids; \
        mv RAxML_bipartitions.tre \$f-prot.tre; rm -f RAxML_bipartitions.idm; done
        cd ../
        mkdir TREE
        cp ali/*.tre TREE/
        echo "GENERA info: Prot ML inference" >> GENERA-Phylogeny.log
        """
    }
    else if (params.mode == 'DNA') {
        println "GENERA info: DNA ML inference"
        """
        mkdir ali
        mv *.ali ali/
        cd ali/
        find *.ali > list
        sed -i -e 's/.ali//g' list
        for f in `cat list`; do ali2phylip.pl \$f.ali --map-ids; done
        for f in `cat list`; do raxmlHPC-PTHREADS-AVX -T $cpu -s \$f.phy \
        -n \$f-RAXML-GTRGAMMA-100xRAPIDBP -m GTRGAMMA -N 100 -f a -x 1975021703574 -p 1975021703574; done
        #format
        for f in `cat list`; do mv \$f.idm RAxML_bipartitions.idm; \
        format-tree.pl RAxML_bipartitions.\$f-RAXML-GTRGAMMA-100xRAPIDBP --map-ids; \
        mv RAxML_bipartitions.tre \$f-DNA.tre; rm -f RAxML_bipartitions.idm; done
        cd ../
        mkdir TREE
        cp ali/*.tre TREE/
        echo "GENERA info: DNA ML inference" >> GENERA-Phylogeny.log
        """
    }
}

//format trees
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'TREE' from mltree_ch1
    file 'IDM' from idm_ch
    file 'ALIGNU' from alignu_ch
    val companion from params.companion
    val version from params.version
    file "GENERA-Phylogeny.log" from log_ch4

    output:
    file 'FORMAT-TREES' into formattreesFINAL_ch
    file 'UNAMBIGUOUS' into unambigousFINAL_ch
    file "GENERA-Phylogeny.log" into logFINAL_ch

    //script
    script:
    println "GENERA info: format-trees"
    """
    mkdir FORMAT-TREES
    mkdir UNAMBIGUOUS
    cp TREE/*.tre FORMAT-TREES/
    cd FORMAT-TREES/
    cp ../IDM .
    find *.tre > list
    sed -i -e 's/.tre//g' list
    for f in `cat list`; do tree2list.pl \$f.tre; sed -i -e 's/ /_/g' \$f.idl; done
    for f in `cat list`; do $companion \$f.idl --mode=IDM; mv idm.temp \$f.idm; done
    for f in `cat list`; do format-tree.pl \$f.tre --map-ids; mv \$f.tre \$f-prot-format.tre; done
    cd ../    
    mv ALIGNU/* UNAMBIGUOUS/
    echo "GENERA info: format-trees" >> GENERA-Phylogeny.log
    echo VERSION: >> GENERA-Phylogeny.log
    echo $version >> GENERA-Phylogeny.log
    """
}
#!/usr/bin/env nextflow

/*
========================================================================================
                         Orthology GENERA
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
    
    Citation:
    Please cite : 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow orthology.sif --protein=prots/
    
    Mandatory arguments:
    --protein               Path to protein sequences files in fasta format.

    Optional arguments:
    --core                  Extract core gene. Default = no.
                            A list of organism is required with this options.
    --list                  Path to a list of organism to consider for core genes.
    --cpu                   Number of cpu to use. Default = 1.                                      

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

//Path to short reads fastq : Mandatory
params.protein = null
if (params.protein == null) {
	exit 1, "Path to protein sequences files in fasta format."
}

//protein
params.core = 'no'

//list
params.list = 'no'

//outdir
params.outdir='GENERA-orthology'

//cpu
params.cpu = '1'

//Path to companion
params.companion = '/opt/COMPANION/Orthology_companion.py'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
protein_ch = Channel.fromPath(params.protein)
list_ch = Channel.fromPath(params.list)


//BMC names
process format {
	//informations

	//input output
    input:
    file '*.faa' from protein_ch

    output:
    file 'FASTA/*abbr.faa' into abbrProteomes_ch1
    file "GENERA-Orthology.log" into log_ch1

    //script
    script:
    println "GENERA info: format names with BMC"
    """
    mkdir FASTA
    cp .faa/*.faa FASTA/
    cd FASTA
    find *.faa > faa.list
    sed -i -e 's/.faa//g' faa.list
    for f in `cat faa.list`; do inst-abbr-ids.pl \$f*.faa --id-regex=:DEF --id-prefix=\$f; done
    cd ../
    echo "GENERA info: format names with BMC" >> GENERA-Orthology.log
    """
}


//Orthology inference
process orthofinder {
	//informations

	//input output
    input:
    file '*-abbr.faa' from abbrProteomes_ch1
    file "GENERA-Orthology.log" from log_ch1
    val cpu from params.cpu

    output:
    file 'OG/OG*.fa' into orthoSeq_ch1
    file 'OG/OG*.fa' into orthoSeq_ch2
    file "GENERA-Orthology.log" into log_ch2

    //script
    script:
    println "GENERA info: run Orthofinder"
    """
    mkdir OF-indir
    mv *.faa OF-indir/
    orthofinder -t $cpu -a $cpu -f OF-indir/
    mkdir OG
    cp OF-indir/OrthoFinder/Results_*/Orthogroup_Sequences/*.fa OG/
    echo "GENERA info: run Orthofinder" >> GENERA-Orthology.log
    """
}

//Orthology inference
process core {
	//informations

	//input output
    input:
    file 'OG/OG*.fa' from orthoSeq_ch1
    file 'list' from list_ch
    file "GENERA-Orthology.log" from log_ch2
    val companion from params.companion

    output:
    file 'CORE/OG*.fa' into core_ch1
    file "GENERA-Orthology.log" into log_ch3


    //script
    script:
    if (params.core == 'yes') {
        if (params.list == 'no') {
            println "GENERA info: organism list must be provided with core genes option"
            """
 
            """
        }
        else {
            println "GENERA info: core gene option activated"
            """
            #mkdir OG
            mkdir CORE
            #mv *.fa OG
            cd OG
            $companion ../list
            cd ../
            mv OG/core-OG.list .
            sed -i -e 's/.fa//g' core-OG.list
            for f in `cat core-OG.list`; do mv OG/\$f.fa CORE/; done
            echo "GENERA info: core gene option activated" >> GENERA-Orthology.log
            echo "Number of Core gene:" >> GENERA-Orthology.log
            wc -l core-OG.list >> GENERA-Orthology.log
            """
        }
    }
    else {
        println "GENERA info: core gene option not activated"
        """
        mkdir CORE
        echo "GENERA info: core gene option not activated" > CORE/OGFALSE.fa
        echo "GENERA info: core gene option not activated" >> GENERA-Orthology.log
        """
    }
}



//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'CORE/OG*.fa' from core_ch1
    file 'OG/OG*.fa' from orthoSeq_ch2
    file "GENERA-Orthology.log" from log_ch3
    val version from params.version


    output:
    file 'coreGenes/OG*.fa' into coreFINAL_ch1
    file 'Orthologous/OG*.fa' into orthoSeqFINAL_ch2
    file "GENERA-Orthology.log" into logFINAL_ch

    //script
    script:
    """
    mkdir coreGenes
    mv CORE/OG*.fa coreGenes/
    mkdir Orthologous
    mv OG/OG*.fa Orthologous/
    echo VERSION: >> GENERA-Orthology.log
    echo $version >> GENERA-Orthology.log
    """
}
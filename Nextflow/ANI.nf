#!/usr/bin/env nextflow

/*
========================================================================================
                        Name of the scripts
========================================================================================
GIT url : 
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

    Version: 1.2.0 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow run ANI.nf --genome=genome --list=list --cpu=20
    
    Mandatory arguments:
    --genome                 Specify the path to genome directory (ext = .fna)
    --list                   Specify the order of organism for the heatmap.

    Optional arguments:
    --mode                   Specifiy the mode, onetomany or manytomany, default == manytomany
    --shortlist              Specify the list of query genomes for onetomany mode
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
params.genome = null
if (params.genome == null) {
	exit 1, "Specify the path to genome directory (ext = .fna)"
}

//Path to OGs : Mandatory
params.list = null
if (params.list == null) {
	exit 1, "Specify the order of organism for the heatmap."
}

//mode
params.mode = 'manytomany'

//short list
params.shortlist = 'none'

//cpu
params.cpu = '1'

//outdir
params.outdir='GENERA_ANI'

//Path to companion
params.companion = '/opt/companion/ANI_companion.py'

//version
params.version = '1.1.0'

/*
CORE PROGRAM
*/

//Load input files
genome_ch = Channel.fromPath(params.genome)
list_ch = Channel.fromPath(params.list)
shortlist_ch1 = Channel.fromPath(params.shortlist)
shortlist_ch2 = Channel.fromPath(params.shortlist)

//BMC names format
process ANI {
	//informations

	//input output
    input:
    file '*' from genome_ch
    file 'shortlist' from shortlist_ch1
    val cpu from params.cpu

    output:
    file 'ANI.txt' into file_ch1
    file 'ANI.txt' into file_ch2

    //script
    script:
    if (params.mode == 'onetomany') {
        println "GENERA info: running fastANI"
        """
        mkdir GEN
        cp genome/* GEN/
        cd GEN/
        find *.fna > list
        cp ../shortlist .
        for f in `cat shortlist`; do fastANI -q \$f.fna --rl list -o ANI-\$f -t $cpu --matrix; done
        rm -f *.matrix
        cd ../
        cat GEN/ANI* > ANI.txt 
        """
    }
    else if (params.mode == 'manytomany') {
        println "GENERA info: running fastANI"
        """
        mkdir GEN
        cp genome/* GEN/
        cd GEN/
        find *.fna > list
        fastANI --ql list --rl list -o ANI -t $cpu --matrix
        cd ../
        cp GEN/ANI ANI.txt 
        """
    }
}

//modellling plots
process heatmap {
	//informations

	//input output
    input:
    file 'ANI.txt' from file_ch1
    file 'list' from list_ch
    file 'shortlist' from shortlist_ch2
    val companion from params.companion

    output:
    file 'ANI.pdf' into pdf_ch
    file 'ANI-close' into close_ch
    file 'ANI-matrix.txt' into matrice_ch

    //script
    script:
    if (params.list == 'none'){
        """
        echo "GENERA info: a list of file should be provided" >> GENERA-METABO.log
        """
    }
    else if (params.mode == 'onetomany') {
        println "GENERA info: heatmap onetomany mode"
        """
        $companion ANI.txt --mode=onetomany
        #Don't run R in onetomany mode, jusr closest list
        echo "onetomany mode" > ANI.pdf
        #close part
        mkdir ANI-close
        mv *closest* ANI-close/
        """        
    }
    else {
        println "GENERA info: heatmap, manytomany mode"
        """
        $companion ANI.txt --mode=manytomany
        #run R
        cp list finallist
        echo "GENOME" > head
        cat head finallist > temp; mv -f temp finallist
        tpage --define file=ANI /opt/code-ANI.tt > code.r
        Rscript code.r
        mv code.r ANI-code.r
        #close part
        mkdir ANI-close
        mv *closest* ANI-close/
        """
    }
}


//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'ANI.txt' from file_ch2
    file 'ANI.pdf' from pdf_ch
    file 'ANI-close' from close_ch
    file 'ANI-matrix.txt' from matrice_ch
    val version from params.version

    output:
    file 'fastANI.txt' into anifinal_ch
    file 'ANI.matrix' into matrixfinal_ch
    file 'heatmap.pdf' into heatmapfinal_ch
    file 'ANI-close-match' into closefinal_ch
    file 'GENERA-ANI.log' into logfinal_ch

    //script
    script:
    """
    cp ANI-matrix.txt ANI.matrix
    cp ANI.txt fastANI.txt
    cp ANI.pdf heatmap.pdf
    cp -r ANI-close/ ANI-close-match
    echo "GENERA info: Running fastANI" >> GENERA-ANI.log
    echo VERSION: >> GENERA-ANI.log
    echo $version >> GENERA-ANI.log
    """
}
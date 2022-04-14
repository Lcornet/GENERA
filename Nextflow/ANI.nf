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

    Version: 1.0.0 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow run ANI.nf --genome=genome --list=list --cpu=20
    
    Mandatory arguments:
    --genome                 Specify the path to genome directory (ext = .fna)

    Optional arguments:
    --list                   Specify the order of organism for the heatmap.
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

//cpu
params.cpu = '1'

//list
params.list = 'none'

//outdir
params.outdir='GENERA_ANI'

//Path to companion
params.companion = '/opt/companion/ANI_companion.py'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
genome_ch = Channel.fromPath(params.genome)
list_ch = Channel.fromPath(params.list)

//BMC names format
process ANI {
	//informations

	//input output
    input:
    file '*' from genome_ch
    val cpu from params.cpu

    output:
    file 'ANI.txt' into file_ch1
    file 'ANI.txt' into file_ch2

    //script
    script:
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

//modellling plots
process heatmap {
	//informations

	//input output
    input:
    file 'ANI.txt' from file_ch1
    file 'list' from list_ch
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
    else {
        println "GENERA info: heatmap"
        """
        $companion ANI.txt
        #run R
        cp list finallist
        echo "GENOME" > head
        cat head finallist > temp; mv -f temp finallist
        tpage --define file=ANI /opt/code-ANI.tt > code.r
        Rscript code.r
        mv code.r ANI-code.r
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
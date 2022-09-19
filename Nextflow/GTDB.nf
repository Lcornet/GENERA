#!/usr/bin/env nextflow

/*
========================================================================================
                        GTDB.nf
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

    nextflow run GTDB.nf --genome=genome --cpu=20
    
    Mandatory arguments:
    --genome                 Specify the path of genome directory (ext, .fna)

    Optional arguments:
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
	exit 1, "Specify the path of genome directory."
}

//cpu
params.cpu = '1'

//GTDBTK_DATA_PATH=/scratch/ulg/GENERA/Databases/GTDB/release207/
params.database='/scratch/ulg/GENERA/Databases/GTDB/release207/'

//outdir
params.outdir='GENERA_GTDB'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
genome_ch = Channel.fromPath(params.genome)


//GTDB
process gtdb {
	//informations

	//input output
    input:
    file '*' from genome_ch
    val cpu from params.cpu
    val database from params.database

    output:
    file 'classify' into gtdb_ch1

    //script
    script:
    """
    mkdir GEN/
    cp genome/*.fna GEN/
    #export GTDBTK_DATA_PATH=/scratch/ulg/GENERA/Databases/GTDB/release207/
    export GTDBTK_DATA_PATH=$database
    gtdbtk classify_wf --genome_dir GEN --out_dir classify_out -x fna --cpus $cpu
    mkdir classify
    cp classify_out/*.tsv classify/
    """
}


//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'classify' from gtdb_ch1
    val version from params.version

    output:
    file 'GTDB-classify' into gtdbfinal_ch
    file 'GENERA-gtdb.log' into logfinal_ch


    //script
    script:
    """
    mkdir GTDB-classify
    cp classify/* GTDB-classify/
    echo "Running GTDB" >> GENERA-gtdb.log
    echo "using release 207 of GTDB" >> GENERA-gtdb.log
    echo VERSION: >> GENERA-gtdb.log
    echo $version >> GENERA-gtdb.log
    """
}
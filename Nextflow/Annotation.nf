#!/usr/bin/env nextflow

/*
========================================================================================
                         Annotation GENERA
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

    nextflow annotation.sif --genome=genome.fasta --mode=prokaryote
    
    Mandatory arguments:
    --genome                    Path to genome in fasta format.
    --mode                      choice between prokayote or eukaryote.

    Optional arguments:
    --organism                  Name of your organism with the format Genus_species or
                                Genus_species_strain (for the augustus gene model name).
                                Mandatory for eukaryotic mode.
    --protein                   Use of protein evidence (boolean switch: 0/1). default= 1.
                                Mandatory for eukaryotic mode.
    --est                       Use of EST/RNA-seq evidence (boolean switch: 0/1). default = 1
                                Mandatory for eukaryotic mode.
    --taxdump                   Path to local mirror of the NCBI Taxonomy database. default = taxdump
                                Mandatory for eukaryotic mode.
    --augustusdb               Path to augustus gene model database. default = augustus-config.
                                Mandatory for eukaryotic mode.
    --augustusgm               Name of an already existing augustus gene model (see the list on the augustus config/species/ folder).
                                Mandatory for eukaryotic mode.
    --protdbs                  Path to the protein databases for the different eukaryotic clades. default = prot_dbs
                                Mandatory for eukaryotic mode. 
    --cpu                       Number of cpus to use. default = 1
    --outdir                    name of the outdir. defaukt = GENERA-annotation

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
params.genome = null
if (params.genome == null) {
	exit 1, "Path to genome in fasta format"
}

params.mode = null
if (params.mode == null) {
	exit 1, "choice between prokayote or eukaryote"
}

//protein
params.protein= '1'

//est
params.est = '1'

//taxdir
params.taxdump = '/scratch/ulg/GENERA/Databases/AMAW/taxdump'

//organism
params.organism = 'none'

//cpu
params.cpu = '1'

//augustus-db
params.augustusdb = '/scratch/ulg/GENERA/Databases/AMAW/augustus-config'

//augustus-gm
params.augustusgm = 'none'

//prot-dbs
params.protdbs = '/scratch/ulg/GENERA/Databases/AMAW/prot_dbs'

//outdir
params.outdir='GENERA-annotation'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
genome_ch = Channel.fromPath(params.genome)
augustusdb_ch = Channel.fromPath(params.augustusdb)
protdbs_ch = Channel.fromPath(params.protdbs)
taxdump_ch = Channel.fromPath(params.taxdump)

//Trimming of short reads
//Download RefSeq metadata, compute taxonomy file and prudce download ftp file
process annotation {
	//informations

	//input output
    input:
    file 'genome.fasta' from genome_ch
    val organism from params.organism
    val protein from params.protein
    val est from params.est
    val taxdump from taxdump_ch
    val cpu from params.cpu
    val augustusdb from augustusdb_ch
    val protdbs from protdbs_ch
    val augustusgm from params.augustusgm

    output:

    //script
    script:
    """
    amaw.pl --genome=genome.fasta --organism=$organism --proteins=$protein --est=$est --taxdir=$taxdump --maker-cpus=$cpu --trinity-cpus=$cpu --rsem-cpus=$cpu --augustus-db=$augustusdb --augustus-gm=$augustusgm --outdir=AMAW --prot-dbs=$protdbs
    """
}

//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:


    output:
    file 'log' into logFINAL_ch

    //script
    script:
    """
    echo TEST > log

    """
}
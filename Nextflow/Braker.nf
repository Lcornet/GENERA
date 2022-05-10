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

    nextflow .... 
    
    Mandatory arguments:
    --genome                 Specify genome
    --currentpath            Specify your current full path (as obtenied by pwd), for TMPDIR

    Optional arguments:
    --brakermode             Specify the mode of braker, with RNAseq + proteins (default = rnaseq) or with protein file only (= prot)
    --SRA                    Specific rnaseq SRA list file, default = none 
    --prot                   Specific which prot file to use, fungi or test, default = fungi 
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
	exit 1, "Specify genome."
}

//Current path : Mandatory
params.currentpath = null
if (params.currentpath == null) {
	exit 1, "Specify the current full path"
}

//brakermode
params.brakermode = 'rnaseq'

//Path to prot
params.prot = 'fungi'

//Path to SRA list
params.SRA = 'testfile'

//cpu
params.cpu = '1'

//outdir
params.outdir='GENERA-braker'

//Path to taxdump
params.config = 'local'

//Not specified by user
//Path to project dir taxdump
confdir = "$workflow.projectDir" + '/Augustus-config'
workingdir = file(confdir)

//Path to companion
params.companion = '/opt/braker_companion.py'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
genome_ch1 = Channel.fromPath(params.genome)
sra_ch1 = Channel.fromPath(params.SRA)
config_ch = Channel.fromPath(params.config)


//Taxonomy, set taxdump if not specifed
process augustusCongig {
	//informations

	//input output
    input:
    val config from config_ch
    
    output:
    file "config_path.txt" into taxdump_path1
    file "config_path.txt" into taxdump_path2
    val confdir into confdir_ch1

    //script
    script:

    taxdir = 'na'

    if (params.config == 'local'){
        println "GENERA-INFO: Config not specified -> project dir"

        if( !workingdir.exists() ) {
            println "GENERA-INFO: Config dir not found in project dir -> Created"
            if( !workingdir.mkdirs() )    {
                exit 1, "Cannot create working config directory"
            }

            confdir = workingdir 

            """
            #mkdir $workingdir
            cp -r /scratch/ulg/GENERA/Databases/BRAKER/Augustus-config-bak/* $workingdir/
            echo $workingdir > config_path.txt
            """
        }
        else {
            println "GENERA-INFO: config dir found in project dir -> Not authorized, erase Augustus-config/"

            confdir = workingdir 

 	        """
            rm -rf $workingdir
            mkdir $workingdir
            cp -r /scratch/ulg/GENERA/Databases/BRAKER/Augustus-config-bak/* $workingdir/
            echo $workingdir > config_path.txt
		    """           

        }
    }

	else{
        println "GENERA-INFO: Config specified"

        confdir = config

		"""
        echo $config > config_path.txt
		"""		
    }
}

//getprot
process getprot {
	//informations

	//input output
    input:

    output:
    file 'proteins.faa' into protcp_ch
    file 'GENERA-braker.log' into log_ch1

    //script
    script:
    if (params.prot == 'test') {
        """
        cp /scratch/ulg/GENERA/Databases/BRAKER/PROTDB/braker-test/proteins.faa proteins.faa
        echo "GENERA info: use test prot file" >> GENERA-braker.log
        """
    }
    else if (params.prot == 'fungi') {
        """
        cp /scratch/ulg/GENERA/Databases/BRAKER/PROTDB/OrthoDB/fungi.faa proteins.faa
        echo "GENERA info: use fungi OrthDB prot file" >> GENERA-braker.log
        """
    }
}

//braker
process abbr {
	//informations

	//input output
    input:
    file 'genome.fna' from genome_ch1
    file 'GENERA-braker.log' from log_ch1

    output:
    file 'genome.fa' into genomeabbr_ch1
    file 'genome.fa' into genomeabbr_ch2
    file 'GENERA-braker.log' into log_ch2

    //script
    script:
    """
    #Genome
    inst-abbr-ids.pl genome.fna --id-regex=:DEF
    mv genome-abbr.fna genome.fa
    sed -i -e 's/>|/>/g' genome.fa
    echo "GENERA info: abbr the genome file" >> GENERA-braker.log
    """
}

//hisat2
process hisat2 {
	//informations

	//input output
    input:
    file 'genome.fna' from genomeabbr_ch1
    file 'SRA' from sra_ch1
    val cpu from params.cpu
    val companion from params.companion
    file 'GENERA-braker.log' from log_ch2

    output:
    file 'BAM' into bam_ch
    file 'bam_list.txt' into bamlist_ch
    file 'GENERA-braker.log' into log_ch3

    //script
    script:
    if (params.brakermode == 'rnaseq') {
        if (params.SRA != 'testfile') {
            """
            for f in `cat SRA`; do fastq-dump --split-files \$f; done
            echo "GENERA info: fastqdump" >> GENERA-braker.log
            hisat2-build genome.fna genome.ht2
            for f in `cat SRA`; do mv \$f*1.fastq \$f-1.fastq; mv \$f*2.fastq \$f-2.fastq; done
            for f in `cat SRA`; do hisat2 -p $cpu -x genome.ht2 -1 \$f-1.fastq -2 \$f-2.fastq -S \$f.sam; done
            for f in `cat SRA`; do samtools view --threads $cpu -b -o \$f.bam \$f.sam; done
            samtools sort -m 10G -o \$f-sorted.bam -T \$f-temp --threads $cpu \$f.bam
            mkdir BAM
            mv *sorted.bam BAM/
            cd BAM
            python3 $companion
            mv bam_list.txt ../
            cd ../
            echo "GENERA info: hisat2 alignment" >> GENERA-braker.log
            """
        }
        else {
            """
            mkdir BAM
            cp /scratch/ulg/GENERA/Databases/BRAKER/PROTDB/braker-test/RNAseq.bam rnaseq.bam
            mv *.bam BAM/
            cd BAM
            python3 $companion
            mv bam_list.txt ../
            cd ../
            echo "GENERA info: bam file provided skipping hisat2" >> GENERA-braker.log
            """
        }
    }
    else {
        """
        mkdir BAM
        echo "GENERA info: no RNAseq file, only prot mode" > BAM/info.bam
        echo "GENERA info: no RNAseq file, only prot mode" > bam_list.txt
        echo "GENERA info: no RNAseq file, only prot mode" >> GENERA-braker.log
        """
    }
}

//braker
process braker {
	//informations
    //publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'proteins.faa' from protcp_ch
    file 'genome.fna' from genomeabbr_ch2
    file 'BAM' from bam_ch
    file 'bam_list.txt' from bamlist_ch
    file "config_path.txt" from taxdump_path1
    file 'GENERA-braker.log' from log_ch3
    val version from params.version
    val path from params.currentpath
    //file 'rnaseq.bam' from rnaseq_ch3

    output:
    file 'braker1_out' into braker_ch
    file 'GENERA-braker.log' into log_ch4


    //script
    script:
    if (params.brakermode == 'rnaseq') {
        """
        export TMPDIR=$path
        cp BAM/* .
        mkdir braker1_out
        braker.pl --genome=genome.fna --bam=\$(<bam_list.txt) --prot_seq=proteins.faa \
        --etpmode --AUGUSTUS_CONFIG_PATH=\$(<config_path.txt) --AUGUSTUS_BIN_PATH=/opt/augustus-3.4.0/bin/ --cores=20  \
        --workingdir=braker1_out 
        #braker.pl --genome=genome.fna --bam=rnaseq.bam --prot_seq=proteins.faa \
        #--etpmode --softmasking --AUGUSTUS_CONFIG_PATH=/scratch/ulg/bioec/lcornet/AMAW/BRAKEN/Augustus-config --AUGUSTUS_BIN_PATH=/opt/augustus-3.4.0/bin/ --cores=20  \
        #--workingdir=braker1_out  
        echo "GENERA info: running braker with RNAseq + proteins" >> GENERA-braker.log
        """
    }
    else {
        """
        export TMPDIR=$path
        mkdir braker1_out
        braker.pl --genome=genome.fna --prot_seq=proteins.faa \
        --AUGUSTUS_CONFIG_PATH=\$(<config_path.txt) --AUGUSTUS_BIN_PATH=/opt/augustus-3.4.0/bin/ --cores=20  \
        --workingdir=braker1_out      
        echo "GENERA info: running braker with proteins only" >> GENERA-braker.log
        """
    }
}

//publications
process results {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'braker1_out' from braker_ch
    file 'GENERA-braker.log' from log_ch4
    val version from params.version

    output:
    file 'BRAKER' into brakerfinal_ch1
    file 'GENERA-braker.log' into logfinal_ch

    //script
    script:
    """
    mv braker1_out BRAKER
    echo "GENERA info: Braker finished" >> GENERA-braker.log
    echo VERSION: >> GENERA-ANI.log
    echo $version >> GENERA-ANI.log
    """
}
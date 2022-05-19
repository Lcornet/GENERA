#!/usr/bin/env nextflow

/*
========================================================================================
                        Metabolic
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

    nextflow run metabolic.nf --infile=infile --mode=modelling --list=list
    
    Mandatory arguments:
    --infile                 Specify the path for proteomes (.faa) genome(s) (.fna) 
    --mode                   Specify the mode, functional (protein) or modelling (genomes)

    Optional arguments:
    --list                   Specify the order of organisnm for the heatmap.
                             Mandatory with modelling mode.    
    --kegg                   Specify the path to kegg dir, automatic build by default.  
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

//Path to genome : Mandatory
params.infile = null
if (params.infile == null) {
	exit 1, "Specify the path with genome(s)"
}

//Path to genome : Mandatory
params.mode = null
if (params.mode == null) {
	exit 1, "Specify the mode, functional or modelling"
}

//list
params.list = 'none'

//cpu
params.cpu = '1'

//outdir
params.outdir='GENERA_metabolic'

//Path to kegg
params.kegg = 'local'

//Not specified by user
//Path to project dir taxdump
keggdir = "$workflow.projectDir" + '/kegg'
workingdir = file(keggdir)

//Path to companion
params.companion = '/opt/companion/metabo_companion.py'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
infile_ch1 = Channel.fromPath(params.infile)
infile_ch2 = Channel.fromPath(params.infile)
list_ch = Channel.fromPath(params.list)
kegg_ch = Channel.fromPath(params.kegg)

//SET KEGG dir with anvio
process keggsetup {
	//informations

	//input output
    input:
    val kegg from kegg_ch
    
    output:
    file "kegg_path.txt" into kegg_path1
    val keggdir into keggdir_ch1

    //script
    script:

    keggdir = 'na'

    if (params.kegg == 'local'){
        println "GENERA-INFO: KEGG dir not specified -> project dir"

        if( !workingdir.exists() ) {
            println "GENERA-INFO: KEGG dir not found in project dir -> Created"

            keggdir = workingdir 

            """
            anvi-setup-kegg-kofams --kegg-data-dir $workingdir
            echo $workingdir > kegg_path.txt
            """
        }
        else {
            println "GENERA-INFO: KEGG dir found in project dir -> Used"

            keggdir = workingdir 

 	        """
            echo $workingdir > kegg_path.txt
		    """           

        }
    }

	else{
        println "GENERA-INFO: Kegg dir specified"

        keggdir = kegg

		"""
        echo $kegg > kegg_path.txt
		"""		
    }
}

//functional
process functional {
	//informations

	//input output
    input:
    file '*' from infile_ch1
    val cpu from params.cpu

    output:
    file 'FUNCT' into mantis_ch
    file "GENERA-METABO.log" into log_ch1


    //script
    script:
    if (params.mode == 'functional') {
        println "GENERA info: functional mode, running mantis"
        """
        #collect protein sequences
        mkdir PROT
        cp infile/*faa PROT
        cd PROT
        #heat *.faa > faalist
        #sed -i -e 's/.faa//g' faalist
        cat *.faa > all.faa
        mantis run -i all.faa -o mantis -ht $cpu -c $cpu
        cd ../
        mkdir FUNCT/
        cp PROT/mantis/* FUNCT/
        echo "GENERA info: functional mode, running mantis" >> GENERA-METABO.log
        """
    }
    else {
        println "GENERA info: functional mode not activated, not running mantis"
        """
        mkdir FUNCT
        echo "GENERA info: functional mode not activated, not running mantis" > FUNCT/info
        echo "GENERA info: functional mode not activated, not running mantis" >> GENERA-METABO.log
        """ 
    }
}

//modelling
process modelling {
	//informations

	//input output
    input:
    file '*' from infile_ch2
    val keggdir from keggdir_ch1
    val cpu from params.cpu
    file "GENERA-METABO.log" from log_ch1

    output:
    file 'MODELLING' into metabo_ch1
    file 'MODELLING' into metabo_ch2
    file "GENERA-METABO.log" into log_ch2

    //script
    script:
    if (params.mode == 'modelling') {
        println "GENERA info: modelling mode, running anvio"
        """
        #collect genomes
        mkdir GEN
        cp infile/*.fna GEN
        #Run anvio
        cd GEN
        find *.fna > fnalist
        sed -i -e 's/.fna//g' fnalist
        for f in `cat fnalist`; do anvi-script-reformat-fasta \$f.fna --simplify-names --seq-type NT \
        --output-file \$f.fa; done
        for f in `cat fnalist`; do anvi-gen-contigs-database --contigs-fasta \$f.fa \
        --project-name \$f --output-db-path \$f-CONTIG.db --num-threads $cpu; done
        echo "GENERA info: modelling mode, running anvio: anvi-run-kegg-kofams" >> GENERA-METABO.log
        for f in `cat fnalist`; do echo \$f; anvi-run-kegg-kofams --kegg-data-dir $keggdir \
        -c \$f-CONTIG.db -T $cpu; done
        for f in `cat fnalist`; do anvi-estimate-metabolism --kegg-data-dir $keggdir \
        -c \$f-CONTIG.db -O \$f-metabo; done
        cd ../
        mkdir MODELLING/
        mv GEN/*metabo_modules.txt MODELLING/
        echo "GENERA info: modelling mode, running anvio" >> GENERA-METABO.log
        """
    }
    else {
        println "GENERA info: modelling mode not activated, not running anvio"
        """
        mkdir MODELLING
        echo "GENERA info: modelling mode not activated, not running anvio" > MODELLING/info
        echo "GENERA info: modelling mode not activated, not running anvio" >> GENERA-METABO.log
        """ 
    }
}

//modellling plots
process modellingplots {
	//informations

	//input output
    input:
    file 'MODELLING' from metabo_ch1
    file 'list' from list_ch
    val companion from params.companion
    file "GENERA-METABO.log" from log_ch2

    output:
    file 'PLOTS' into plots_ch1
    file "GENERA-METABO.log" into log_ch3

    //script
    script:
    if (params.mode == 'modelling') {
        if (params.list == 'none'){
            """
            echo "GENERA info: a list of file should be provided with medlloing mode" >> GENERA-METABO.log
            """
        }
        else {
            println "GENERA info: modelling mode, making plots"
            """
            cp MODELLING/* .
            $companion list --mode=no-zero
            find *infile.txt > infile.list
            sed -i -e 's/-infile.txt//g' infile.list
            #run R
            cp list finallist
            echo "GENOME" > head
            cat head finallist > temp; mv -f temp finallist
            for f in `cat infile.list`; do tpage --define file=\$f /opt/code.tt > code.r; \
            Rscript code.r; mv code.r \$f-code.r ; done
            #merge pdf
            pdftk *.pdf cat output merge.pdf
            mkdir PLOTS
            mv *.pdf PLOTS/
            echo "GENERA info: modelling mode, making plots" >> GENERA-METABO.log
            """
        }
    }
    else {
        println "GENERA info: modelling mode not activated, no plots"
        """
        mkdir PLOTS
        echo "GENERA info: modelling mode not activated, no plots" > PLOTS/info
        echo "GENERA info: modelling mode not activated, no plots" >> GENERA-METABO.log
        """ 
    }
}

//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'MODELLING' from metabo_ch2
    file 'FUNCT' from mantis_ch
    file 'PLOTS' from plots_ch1
    val version from params.version
    file "GENERA-METABO.log" from log_ch3

    output:
    file 'ANVIO' into anviofinal_ch1
    file 'MANTIS' into mantisfinal_ch1
    file 'PDF' into plotsfinal_ch1
    file 'GENERA-METABO.log' into finallog_ch

    //script
    script:
    """
    #functional part
    mkdir MANTIS
    mv FUNCT/* MANTIS
    #modelling part
    mkdir ANVIO
    mv MODELLING/* ANVIO
    mkdir PDF
    mv PLOTS/* PDF/
    echo VERSION: >> GENERA-METABO.log
    echo $version >> GENERA-METABO.log
    """
}
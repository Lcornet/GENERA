#!/usr/bin/env nextflow

/*
========================================================================================
                         Orthology
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

    Version: 2.0.6 
    
    Citation:
    Please cite : 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow run Orthology.nf --mode=inference --infiles=infiles --type=nucleotide \
    --core=yes --corelist=corelist --specific=yes --specificlist=specificlist --anvio=yes 
    
    Mandatory arguments:
    --mode                   Specify mode, OG or inference mode.  
    --infiles                Path to sequences files: '.fna' for nucleotide, or '.faa' for protein or '.fa' for OG.
    --type                   Specify infile type: nucleotide (only for prokaryote), protein. 

    Optional arguments:
    --anvio                  Activate anvio pangenomic for orthology, only available for bacteria, yes or no. Default = no.
                             This option will also determine the origin of OGs (Orthofinder or anvio) when the mode is set as OG.
    --core                   Extract core gene. Default = no.
                             A list of organism is required with this options.
    --corelist               Path to a list of organism to consider for core genes.
    --coreunwanted           Specify the maximal number of organism outside of corelist to consider an OGs as core; default = 0
    --coreminfunctionalindex Specify the minimal anvio functional index to consider an OGs as core; default = 0.8
    --coremingeometricindex  Specify the minimal anvio geometric index to consider an OGs as core; default = 0.8
    --specific               Extract specific gene. Default = no.
                             A list of organism is required with this options.
    --specificlist           Path to a list of organism to consider for specific genes.
    --progigalProcedure      Select procedure of prodigal, meta or single, default = single
    --taxdump 		         Path to taxdump directory - automatic setup by default
    --cpu                    Number of cpu to use. Default = 1.  
    --COG                    Path to COG as setup by anvio, default = none                                  

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

//Path to list
//params.corelist = null
//if (params.corelist == null) {
//	exit 1, "Path to a list of organism to consider for specific genes."
//}

//mode
params.mode = null
if (params.mode == null) {
	exit 1, "Specify mode, OG or inference mode."
}

//infiles
params.infiles = null
if (params.infiles == null) {
	exit 1, "Path to sequences files: .faa for protein or .fa for OGs."
}

//infiles
params.type = null
if (params.type == null) {
	exit 1, "Specify infile type: nucleotide (only for prokaryote), protein or OG."
}

//anvio
params.anvio = 'no'

//core
params.core = 'no'

//Path to core list
params.corelist = 'none'

//core unwanted
params.coreunwanted = '0'

//core coreminfunctionalindex
params.coreminfunctionalindex = '0.8'

//core coremaxgeometricindex
params.coremingeometricindex = '0.8'

//specific
params.specific = 'no'

//Path to list
params.specificlist = 'none'

//progigalProcedure
params.progigalProcedure = 'single'

//Path to COG
params.COG = '/scratch/ulg/GENERA/Databases/ANVIO/COG/'

//outdir
params.outdir='GENERA-Orthology'

//cpu
params.cpu = '1'

//Path to taxdump
params.taxdump = 'local'

//Not specified by user
//Path to project dir taxdump
taxdir = "$workflow.projectDir" + '/taxdump'
workingdir = file(taxdir)

//Path to companion
params.companion = '/opt/COMPANION/Orthology_companion.py'
params.anviopantoogs = '/opt/anvio_pan-to-OGs.py'
params.anvioOGsfiltration = '/opt/anvio_OGs-filtration.py'
params.anviotoBMC = '/opt/anvio-to-BMC.py'
params.confirmog = '/opt/confirm-OG.py'
params.ftcompanion = '/opt/companion/OGsEnrichment_companion.py'
params.changespadesids = '/opt/change-spades-IDs.py'

//version
params.version = '2.0.6'

/*
CORE PROGRAM
*/

//Load input files
taxdump_ch = Channel.fromPath(params.taxdump)
infiles_ch = Channel.fromPath(params.infiles)
corelist_ch = Channel.fromPath(params.corelist)
specificlist_ch = Channel.fromPath(params.specificlist)
specificlist_ch2 = Channel.fromPath(params.specificlist)
specificlist_ch3 = Channel.fromPath(params.specificlist)

//Taxonomy, set taxdump if not specifed
process Taxonomy {
	//informations

	//input output
    input:
    val taxdump from taxdump_ch
    
    output:
    file "taxdump_path.txt" into taxdump_path1
    val taxdir into taxdir_ch1

    //script
    script:

    taxdir = 'na'

    if (params.taxdump == 'local'){
        println "GENERA-INFO: Taxdump not specified -> project dir"

        if( !workingdir.exists() ) {
            println "GENERA-INFO: Taxdump dir not found in project dir -> Created"
            if( !workingdir.mkdirs() )    {
                exit 1, "Cannot create working directory"
            }

            taxdir = workingdir 

            """
            setup-taxdir.pl --taxdir=$workingdir
            echo $workingdir > taxdump_path.txt
            """
        }
        else {
            println "GENERA-INFO: Taxdump dir found in project dir -> Used"

            taxdir = workingdir 

 	        """
            echo $workingdir > taxdump_path.txt
		    """           

        }
    }

	else{
        println "GENERA-INFO: Taxdump specified"

        taxdir = taxdump

		"""
        echo $taxdump > taxdump_path.txt
		"""		
    }
}

//BMC names
process format {
	//informations

	//input output
    input:
    file '*' from infiles_ch
    val changespadesids from params.changespadesids

    output:
    file 'FORMAT' into seqs_ch1
    file 'FORMAT' into seqs_ch2
    file "GENERA-SGC.log" into log_ch1

    //script
    script:
    if (params.mode == 'inference') {
        println "GENERA info: inference mode, format names with BMC"
        if (params.type == 'protein') {
            println "GENERA info: inference mode, format protein names with BMC"
            """
            mkdir FASTA
            cp infiles/*.faa FASTA/
            cd FASTA
            find *.faa > faa.list
            sed -i -e 's/.faa//g' faa.list
            for f in `cat faa.list`; do inst-abbr-ids.pl \$f*.faa --id-regex=:DEF --id-prefix=\$f; done
            cd ../
            mkdir FORMAT
            mv FASTA/*abbr.faa FORMAT/
            #as for anvio, delete the -abbr
            cd FORMAT/
            find *.faa > abbr.list
            sed -i -e 's/-abbr.faa//g' abbr.list
            for f in `cat abbr.list`; do mv \$f-abbr.faa \$f.faa; done
            cd ../
            echo "GENERA info: inference mode, format nucleotide names with BMC" >> GENERA-SGC.log
            """
        }
        else if (params.type == 'nucleotide') {
            println "GENERA info: inference mode, format nucleotide names with BMC"
            """
            mkdir CHANGE
            cp infiles/*.fna CHANGE/
            cd CHANGE/
            find *.fna > fna.list
            sed -i -e 's/.fna//g' fna.list
            $changespadesids
            for f in `cat fna.list`; do mv \$f-spades.fna \$f.fna; done
            cd ../
            mkdir FASTA
            mv CHANGE/*.fna FASTA/
            cd FASTA
            find *.fna > fna.list
            sed -i -e 's/.fna//g' fna.list
            for f in `cat fna.list`; do inst-abbr-ids.pl \$f*.fna --id-regex=:DEF --id-prefix=\$f; done
            cd ../
            mkdir FORMAT
            mv FASTA/*abbr.fna FORMAT/
            # For anvio run, reomve the -abbr
            cd FORMAT/
            find *.fna > abbr.list
            sed -i -e 's/-abbr.fna//g' abbr.list
            for f in `cat abbr.list`; do mv \$f-abbr.fna \$f.fna; done
            cd ../
            echo "GENERA info: inference mode, format nucleotide names with BMC" >> GENERA-SGC.log
            """
        }
    }
    else if (params.mode == 'OG') {
        println "GENERA info: OG mode, no formating"
        """
        mkdir FORMAT/
        cp infiles/* FORMAT/
        echo "GENERA info: OG mode, no formating" >> GENERA-SGC.log
        """
    }
}

//prodigal
process prodigal {
	//informations

	//input output
    input:
    file 'FORMAT' from seqs_ch1
    file "GENERA-SGC.log" from log_ch1
    val progigalProcedure from params.progigalProcedure

    output:
    file 'PSEQS' into prot_ch1
    file 'PSEQS' into prot_ch2
    file "GENERA-SGC.log" into log_ch2

    //script
    script:
    if (params.mode == 'inference') {
        if (params.type == 'protein') {
            println "GENERA info: protein mode, will not run prodigal"
            """
            mkdir PSEQS
            mv FORMAT/*.faa PSEQS/
            echo "GENERA info: protein mode, will not run prodigal" >> GENERA-SGC.log
            """
        }
        else if (params.type == 'nucleotide') {
            println "GENERA info: nucleotide mode, run prodigal"
            """
            mkdir PSEQS
            cd FORMAT/
            find *.fna > list
            sed -i -e 's/.fna//g' list
            for f in `cat list`; do prodigal -i \$f.fna -o \$f.genes -d \$f.genes.fna -a \$f.faa -p $progigalProcedure; done
            rm -rf *genes*
            #abbr the file and delete the -abbr in the file nams
            for f in `cat list`; do inst-abbr-ids.pl \$f.faa --id-regex=:DEF; done
            for f in `cat list`; do mv \$f-abbr.faa \$f.faa; done
            sed -i -e 's/>|/>/g' *.faa
            cd ../
            mv FORMAT/*.faa PSEQS/
            echo "GENERA info: nucleotide mode, run prodigal" >> GENERA-SGC.log
            """
        }
    }
    else if (params.mode == 'OG') {
        println "GENERA info: OG mode, will not run prodigal"
        """
        mkdir PSEQS
        mv FORMAT/* PSEQS/
        echo "GENERA info: OG mode, no formating" >> GENERA-SGC.log
        """
    }
}

//Orthology inference
process orthofinder {
	//informations

	//input output
    input:
    file 'PSEQS' from prot_ch1
    val cpu from params.cpu
    file "GENERA-SGC.log" from log_ch2

    output:
    file 'OG' into orthoSeq_ch1
    file "GENERA-SGC.log" into log_ch3

    //script
    script:
    if (params.mode == 'inference') {
        if (params.anvio == 'yes') {
            println "GENERA info: anvio activated, will not run OF"
            """
            mkdir OG
            echo "GENERA info: anvio activated, will not run OF" > OG/FALSE.fa
            echo "GENERA info: anvio activated, will not run OF" >> GENERA-SGC.log
            """           
        }
        else {
            println "GENERA info: run Orthofinder"
            """
            mkdir OF-indir
            cp PSEQS/*.faa OF-indir/
            orthofinder -t $cpu -a $cpu -f OF-indir/
            mkdir OG
            cp OF-indir/OrthoFinder/Results_*/Orthogroup_Sequences/*.fa OG/
            echo "GENERA info: run Orthofinder" >> GENERA-SGC.log
            """
        }
    }
    else if (params.mode == 'OG') {
        println "GENERA info: OG mode, will not run Orthology"
        """
        mkdir OG
        cp PSEQS/* OG/
        echo "GENERA info: OG mode, will not run Orthofinder" >> GENERA-SGC.log
        """
    }
}

//anvio
process anvio {
	//informations

	//input output
    input:
    file 'FORMAT' from seqs_ch2
    val cpu from params.cpu
    val cog from params.COG
    file "GENERA-SGC.log" from log_ch3

    output:
    file 'ANVIO' into anviopan_ch1
    file 'ANVIO' into anviopan_ch2
    file 'ANVIO' into anviopan_ch3
    file "GENERA-SGC.log" into log_ch4

    //script
    script:
    if (params.mode == 'inference') {
        if (params.anvio == 'yes') {
            if (params.type == 'nucleotide') {
                println "GENERA info: run anvio"
                """
                mkdir anvio
                mv FORMAT/*.fna anvio/
                cd anvio/
                find *.fna > list
                sed -i -e 's/.fna//g' list
                #for f in `cat list`; do sed -i -e 's/|/-/g' \$f.fna; done
                for f in `cat list`; do anvi-script-reformat-fasta \$f.fna -o \$f-fixed.fna -l 0 --simplify-names --seq-type NT; done
                for f in `cat list`; do anvi-gen-contigs-database -f \$f-fixed.fna -o \$f.db -n \$f --num-threads $cpu; done
                for f in `cat list`; do anvi-run-ncbi-cogs -c \$f.db --cog-data-dir $cog -T $cpu; done
                find *.db > db.list
                cut -f1 -d'.' db.list > name.list
                echo 'name' > temp
                cat temp name.list > f1
                echo 'contigs_db_path' > temp
                cat temp db.list > f2
                paste f1 f2 > external-genomes.txt
                for f in `cat list`; do anvi-gen-genomes-storage -e external-genomes.txt -o GENOMES.db; done
                anvi-pan-genome -g GENOMES.db --project-name 'PAN' --output-dir PAN-anvio --num-threads $cpu --mcl-inflation 10 --min-occurrence 2
                anvi-get-sequences-for-gene-clusters -g GENOMES.db -p PAN-anvio/PAN-PAN.db -o PAN.fasta
                anvi-script-add-default-collection -p PAN-anvio/PAN-PAN.db -C COL
                anvi-summarize -g GENOMES.db -p PAN-anvio/PAN-PAN.db -o SUMMARY -C COL
                anvi-compute-gene-cluster-homogeneity -g GENOMES.db -p PAN-anvio/PAN-PAN.db -o homogeneity.txt 
                cd ../
                mkdir ANVIO
                mv anvio/SUMMARY/PAN_gene_clusters_summary.txt.gz ANVIO/
                mv anvio/homogeneity.txt ANVIO/
                mv anvio/PAN.fasta ANVIO/
                mv anvio/external-genomes.txt ANVIO/
                echo "GENERA info: run anvio" >> GENERA-SGC.log
                """  
            }
            else {
                println "GENERA info: anvio is only accessibe for nucleotide: will die"
                """
                echo "GENERA info: anvio is only accessibe for nucleotide: will die" >> GENERA-SGC.log
                """ 
            }         
        }
        else {
            println "GENERA info: anvio not activated, OF already run"
            """
            mkdir ANVIO
            echo "GENERA info: anvio not activated, OF already run" >> ANVIO/info.file
            echo "GENERA info: anvio not activated, OF already run" >> GENERA-SGC.log
            """
        }
    }
    else if (params.mode == 'OG') {
        println "GENERA info: OG mode, will not run Orthology inference"
        """
        mkdir ANVIO
        echo "GENERA info: OG mode, will not run Orthology inference" >> ANVIO/info.file
        echo "GENERA info: OG mode, will not run Orthology inference" >> GENERA-SGC.log
        """
    }
}


//Specific gene content
process formatOG {
	//informations

	//input output
    input:
    file 'ANVIO' from anviopan_ch1
    file 'OG' from orthoSeq_ch1
    val companion from params.companion
    val anviopantoogs from params.anviopantoogs 
    val anviotoBMC from params.anviotoBMC
    file "GENERA-SGC.log" from log_ch4

    output:
    file 'ORTHO' into ortho_ch1
    file 'ORTHO' into ortho_ch2
    file "GENERA-SGC.log" into log_ch5


    //script
    script:
    if (params.mode == 'inference') {
        if (params.anvio == 'yes') {
            println "GENERA info: anvio files"
            """
            mkdir ORTHO
            cd ANVIO/
            #anvio OGs
            $anviopantoogs PAN.fasta
            #BMC OGs
            cut -f2 external-genomes.txt > f1
            sed -i -e 's/contigs_db_path/#/g' f1
            sed -i -e 's/.db//g' f1
            cut -f1 external-genomes.txt > f2
            sed -i -e 's/name/#/g' f2
            paste f1 f2 > IDM
            $anviotoBMC --idm=IDM
            cd ../
            mv ANVIO/*OG.fasta ORTHO/
            echo "GENERA info: anvio files" >> GENERA-SGC.log
            """         
        }
        else {
            println "GENERA info: anvio not activated, Orthofinder files"
            """
            mkdir ORTHO
            mv OG/*fa ORTHO/
            echo "GENERA info: anvio not activated, Orthofinder files" >> GENERA-SGC.log
            """
        }
    }
    else if (params.mode == 'OG') {
        println "GENERA info: OG mode, no anvio or orthofinder, copy of OGs"
        """
        mkdir ORTHO/
        mv OG/* ORTHO/
        echo "GENERA info: OG mode, no anvio or orthofinder, copy of OGs" >> GENERA-SGC.log
        """
    }
}

//Core gene content
process core {
	//informations

	//input output
    input:
    file 'ORTHO' from ortho_ch1
    file 'ANVIO' from anviopan_ch2
    file 'list' from corelist_ch
    val companion from params.companion
    val anvioOGsfiltration from params.anvioOGsfiltration
    val coreunwanted from params.coreunwanted 
    val coreminfunctionalindex from params.coreminfunctionalindex
    val coremingeometricindex from params.coremingeometricindex
    file "GENERA-SGC.log" from log_ch5

    output:
    file 'core-OG.list' into coreOGlist_ch
    file 'CORE' into core_ch
    file 'Redo' into redo1_ch
    file "GENERA-SGC.log" into log_ch6


    //script
    script:
    if (params.core == 'yes') {
        if (params.anvio == 'yes') {
            if (params.mode == 'inference') {
                println "GENERA info: core gene from anvio files"
                """
                mkdir OG-anvio
                mkdir OG-BMC
                mkdir Redo
                cp ORTHO/*BMC_OG.fasta OG-BMC/
                cp ORTHO/*OG.fasta OG-anvio/
                cp ANVIO/PAN_gene_clusters_summary.txt.gz OG-anvio/
                cp ANVIO/homogeneity.txt OG-anvio/
                cp OG-anvio/* Redo/
                #Anvio don't keep word after points, do the same thing in the list
                cut -f1 -d"." list > corelist
                cd OG-anvio/
                rm -rf *BMC*
                gunzip PAN_gene_clusters_summary.txt.gz
                $anvioOGsfiltration ../corelist --pfilter=yes --fraction=1 --unwanted=$coreunwanted --cfilter=yes --maxcopy=1 \
                --hfilter=yes --hindex=homogeneity.txt --maxfunctionalindex=$coreminfunctionalindex --maxgeometricindex=$coremingeometricindex \
                --cogtable=PAN_gene_clusters_summary.txt
                cd ../
                mkdir CORE
                mv OG-anvio/filtered-OG.list .
                grep -v "#" filtered-OG.list | grep -v 'OG-name' | cut -f1 > core-OG.list
                for f in `cat core-OG.list`; do cp OG-BMC/\$f*.fasta CORE/; done
                echo "Number of Core gene:" >> GENERA-SGC.log
                wc -l core-OG.list >> GENERA-SGC.log
                echo "GENERA info: core gene from anvio files" >> GENERA-SGC.log
                """ 
            }
            else if (params.mode == 'OG') {
                println "GENERA info: core gene from pre-computed anvio files"
                """
                mkdir OG-anvio
                mkdir OG-BMC
                mkdir Redo
                cp ORTHO/*BMC_OG.fasta OG-BMC/
                cp ORTHO/*OG.fasta OG-anvio/
                cp ORTHO/PAN_gene_clusters_summary.txt.gz OG-anvio/
                cp ORTHO/homogeneity.txt OG-anvio/
                cp OG-anvio/* Redo/
                #Anvio don't keep word after points, do the same thing in the list
                cut -f1 -d"." list > corelist
                cd OG-anvio/
                rm -rf *BMC*
                gunzip PAN_gene_clusters_summary.txt.gz
                $anvioOGsfiltration ../corelist --pfilter=yes --fraction=1 --unwanted=$coreunwanted --cfilter=yes --maxcopy=1 \
                --hfilter=yes --hindex=homogeneity.txt --maxfunctionalindex=$coreminfunctionalindex --maxgeometricindex=$coremingeometricindex \
                --cogtable=PAN_gene_clusters_summary.txt
                cd ../
                mkdir CORE
                mv OG-anvio/filtered-OG.list .
                grep -v "#" filtered-OG.list | grep -v 'OG-name' | cut -f1 > core-OG.list
                for f in `cat core-OG.list`; do cp OG-BMC/\$f*.fasta CORE/; done
                echo "Number of Core gene:" >> GENERA-SGC.log
                wc -l core-OG.list >> GENERA-SGC.log
                echo "GENERA info: core gene from anvio files" >> GENERA-SGC.log
                """   
            }      
        }
        else {
            if (params.mode == 'inference') {
                println "GENERA info: core gene from Orthofinder files"
                """
                mkdir CORE
                mkdir OG-BMC/
                mkdir Redo/
                cp ORTHO/*.fa OG-BMC/
                cp OG-BMC/* Redo/
                cd OG-BMC/
                $companion ../list --unwanted=$coreunwanted
                cd ../
                mv OG-BMC/core-OG.list .
                sed -i -e 's/.fa//g' core-OG.list
                for f in `cat core-OG.list`; do cp OG-BMC/\$f.fa CORE/; done
                echo "Number of Core gene:" >> GENERA-SGC.log
                wc -l core-OG.list >> GENERA-SGC.log
                echo "GENERA info: core gene from Orthofinder files" >> GENERA-SGC.log
                """
            }
            else if (params.mode == 'OG') {
                println "GENERA info: core gene from pre-computed Orthofinder files"
                """
                mkdir CORE
                mkdir OG-BMC/
                mkdir Redo/
                cp ORTHO/*.fa OG-BMC/
                cp OG-BMC/* Redo/
                cd OG-BMC/
                $companion ../list --unwanted=$coreunwanted
                cd ../
                mv OG-BMC/core-OG.list .
                sed -i -e 's/.fa//g' core-OG.list
                for f in `cat core-OG.list`; do cp OG-BMC/\$f.fa CORE/; done
                echo "Number of Core gene:" >> GENERA-SGC.log
                wc -l core-OG.list >> GENERA-SGC.log
                echo "GENERA info: core gene from pre-computed Orthofinder files" >> GENERA-SGC.log
                """
            }
        }
    }
    else {
        println "GENERA info: core gene option not activated"
        """
        mkdir CORE
        mkdir REDO
        echo "GENERA info: core gene option not activate" > core-OG.list
        echo "GENERA info: core gene option not activated" >> CORE/info
        echo "GENERA info: core gene option not activated: REDO option available only after core genes" >> REDO/info
        echo "GENERA info: core gene option not activated" >> GENERA-SGC.log  
        """
    }
}

//Specific gene content
process specific {
	//informations

	//input output
    input:
    file 'ORTHO' from ortho_ch2
    file 'ANVIO' from anviopan_ch3
    file 'list' from specificlist_ch
    val companion from params.companion
    val anvioOGsfiltration from params.anvioOGsfiltration
    file "GENERA-SGC.log" from log_ch6

    output:
    file 'CANDIDATEspec' into candidate_ch
    file 'ANVIO-files' into anviofiles_ch
    file 'OG-BMC' into ogbmc_ch
    file "GENERA-SGC.log" into log_ch7


    //script
    script:
    if (params.specific == 'yes') {
        if (params.anvio == 'yes') {
            if (params.mode == 'inference') {
                println "GENERA info: specific gene from anvio files"
                """
                mkdir OG-anvio
                mkdir OG-BMC
                mv ORTHO/*BMC_OG.fasta OG-BMC/
                mv ORTHO/*OG.fasta OG-anvio/
                mv ANVIO/PAN_gene_clusters_summary.txt.gz OG-anvio/
                mv ANVIO/homogeneity.txt OG-anvio/
                #Anvio don't keep word after points, do the same thing in the list
                cut -f1 -d"." list > specificlist
                cd OG-anvio/
                gunzip PAN_gene_clusters_summary.txt.gz
                $anvioOGsfiltration ../specificlist --pfilter=yes --fraction=0.6 --unwanted=0 --cfilter=yes --maxcopy=1 \
                --hfilter=yes --hindex=homogeneity.txt --maxfunctionalindex=0.8 --maxgeometricindex=0.8 \
                --cogtable=PAN_gene_clusters_summary.txt
                cd ../
                mkdir CANDIDATEspec
                mv OG-anvio/filtered-OG.list .
                grep -v "#" filtered-OG.list | grep -v 'OG-name' | cut -f1 > specific-OG.list
                for f in `cat specific-OG.list`; do cp OG-BMC/\$f*.fasta CANDIDATEspec/; done
                #get files
                mkdir ANVIO-files
                cp OG-anvio/homogeneity.txt ANVIO-files/
                cp OG-anvio/PAN_gene_clusters_summary.txt ANVIO-files/
                cp filtered-OG.list ANVIO-files/
                echo "GENERA info: specific gene from anvio files" >> GENERA-SGC.log
                """ 
            }
            else if (params.mode == 'OG') {
                println "GENERA info: specific gene from pre-computed anvio files"
                """
                mkdir OG-anvio
                mkdir OG-BMC
                mv ORTHO/*BMC_OG.fasta OG-BMC/
                mv ORTHO/*OG.fasta OG-anvio/
                mv ORTHO/PAN_gene_clusters_summary.txt OG-anvio/
                mv ORTHO/homogeneity.txt OG-anvio/
                #Anvio don't keep word after points, do the same thing in the list
                cut -f1 -d"." list > specificlist
                cd OG-anvio/
                $anvioOGsfiltration ../specificlist --pfilter=yes --fraction=0.6 --unwanted=0 --cfilter=yes --maxcopy=1 \
                --hfilter=yes --hindex=homogeneity.txt --maxfunctionalindex=0.8 --maxgeometricindex=0.8 \
                --cogtable=PAN_gene_clusters_summary.txt
                cd ../
                mkdir CANDIDATEspec
                mv OG-anvio/filtered-OG.list .
                grep -v "#" filtered-OG.list | grep -v 'OG-name' | cut -f1 > specific-OG.list
                for f in `cat specific-OG.list`; do cp OG-BMC/\$f*.fasta CANDIDATEspec/; done
                #get files
                mkdir ANVIO-files
                echo "GENERA info: specific gene from anvio files >> ANVIO-files/info
                echo "GENERA info: specific gene from anvio files" >> GENERA-SGC.log
                """   
            }      
        }
        else {
            if (params.mode == 'inference') {
                println "GENERA info: specific gene from Orthofinder files"
                """
                mkdir CANDIDATEspec
                mkdir OG-BMC/
                mv ORTHO/*.fa OG-BMC/
                cd OG-BMC/
                $companion ../list --presence=60 --unwanted=0
                cd ../
                mv OG-BMC/core-OG.list specific-OG.list
                sed -i -e 's/.fa//g' specific-OG.list
                for f in `cat specific-OG.list`; do cp OG-BMC/\$f.fa CANDIDATEspec/\$f.fasta; done
                echo "GENERA info: specific gene from Orthofinder files" > filtered-OG.list
                #get files
                mkdir ANVIO-files
                echo "GENERA info: specific gene from Orthofinder files" >> ANVIO-files/info
                echo "GENERA info: specific gene from Orthofinder files" >> GENERA-SGC.log
                """
            }
            else if (params.mode == 'OG') {
                println "GENERA info: specific gene from pre-computed Orthofinder files"
                """
                mkdir CANDIDATEspec
                mkdir OG-BMC/
                mv ORTHO/*.fa OG-BMC/
                cd OG-BMC/
                $companion ../list --presence=60 --unwanted=0
                cd ../
                mv OG-BMC//core-OG.list specific-OG.list
                sed -i -e 's/.fa//g' specific-OG.list
                for f in `cat specific-OG.list`; do cp OG-BMC/\$f.fa CANDIDATEspec/\$f.fasta; done
                echo "GENERA info: specific gene from Orthofinder files" > filtered-OG.list
                #get files
                mkdir ANVIO-files
                echo "GENERA info: specific gene from pre-computed Orthofinder files >> ANVIO-files/info
                echo "GENERA info: specific gene from pre-computed Orthofinder files" >> GENERA-SGC.log
                """
            }
        }
    }
    else {
        if (params.anvio == 'yes') {
            println "GENERA info: core gene option not activated"
            """
            mkdir CANDIDATEspec
            mkdir ANVIO-files
            mkdir OG-BMC
            mv ORTHO/*BMC_OG.fasta OG-BMC/
            echo "GENERA info: specific gene option not activated" > filtered-OG.list
            echo "GENERA info: specific gene option not activated" >> CANDIDATEspec/info
            echo "GENERA info: specific gene option not activated" >> ANVIO-files/info
            echo "GENERA info: specific gene option not activated" >> OG-BMC/info
            echo "GENERA info: specific gene option not activated" >> GENERA-SGC.log  
            """
        }
        else {
            println "GENERA info: core gene option not activated"
            """
            mkdir CANDIDATEspec
            mkdir ANVIO-files
            mkdir OG-BMC
            mv ORTHO/*.fa OG-BMC/
            echo "GENERA info: specific gene option not activated" > filtered-OG.list
            echo "GENERA info: specific gene option not activated" >> CANDIDATEspec/info
            echo "GENERA info: specific gene option not activated" >> ANVIO-files/info
            echo "GENERA info: specific gene option not activated" >> OG-BMC/info
            echo "GENERA info: specific gene option not activated" >> GENERA-SGC.log  
            """            
        }
    }
}

//Orthologous enrichment part
process enrichment {
	//informations

	//input output
    input:
    file 'PSEQS' from prot_ch2
    file 'CANDIDATEspec' from candidate_ch
    file 'list' from specificlist_ch2
    val ftcompanion from params.ftcompanion
    val cpu from params.cpu
    val taxdir from taxdir_ch1
    file "GENERA-SGC.log" from log_ch7

    output:
    file 'ENRICHED' into enriched_ch
    file "GENERA-SGC.log" into log_ch8

    //script
    script:
    if (params.specific == 'yes') {
        println "GENERA info: running specific enrichment"
        """
        echo "GENERA info: running specific enrichment" >> GENERA-SGC.log

        #Reference organism part
        mkdir ref-banks
        for f in `cat list`; do cp PSEQS/\$f*.faa ref-banks/\$f.faa; done
        cd ref-banks/
        find *.faa > list
        sed -i -e 's/.faa//g' list
        for f in `cat list`; do makeblastdb -in \$f.faa -dbtype prot -parse_seqids -out \$f; done
        paste list list > ref-bank-mapper.idm
        cd ../
        echo "GENERA info: 42, representative bank done" >> GENERA-SGC.log

        #Define Queries
        #cp ref-banks/list .
        $ftcompanion list --mode=queries
        echo "GENERA info: 42, queries done" >> GENERA-SGC.log

        #Part for org to add
        mkdir org-to-add
        cd PSEQS/
        find *.faa > abbr.list
        sed -i -e 's/.faa//g' abbr.list
        cd ../
        cp PSEQS/abbr.list .
        for f in `cat abbr.list`; do cp PSEQS/\$f.faa org-to-add/\$f.faa; done
        cd org-to-add/
        find *.faa > list
        sed -i -e 's/.faa//g' list
        for f in `cat list`; do makeblastdb -in \$f.faa -dbtype prot -parse_seqids -out \$f; done
        paste list list > bank-mapper.idm
        cd ..
        echo "GENERA info: 42, org bank done" >> GENERA-SGC.log

        #Part OGs
        mkdir ORTHO
        cp CANDIDATEspec/*.fasta ORTHO/
        cd ORTHO/
        find *.fasta > list
        sed -i -e 's/.fasta//g' list
        for f in `cat list`; do fasta2ali.pl \$f.fasta; done
        cd ../
        echo "GENERA info: 42, OGs done" >> GENERA-SGC.log

        #Generate yaml
        yaml-generator-42.pl --run_mode=phylogenomic --out_suffix=-GENERA --queries queries.idl --evalue=1e-05 --homologues_seg=yes \
        --max_target_seqs=10000 --templates_seg=no --bank_dir org-to-add --bank_suffix=.psq --bank_mapper org-to-add/bank-mapper.idm --ref_brh=on \
        --ref_bank_dir ref-banks --ref_bank_suffix=.psq --ref_bank_mapper ref-banks/ref-bank-mapper.idm --ref_org_mul=0.3 --ref_score_mul=0.99 \
        --trim_homologues=on --ali_keep_lengthened_seqs=keep --aligner_mode=blast --tax_reports=off --tax_dir $taxdir \
        --megan_like --tol_check=off
        echo "GENERA info: 42, Yaml done" >> GENERA-SGC.log

        #run forty-two
        forty-two.pl ORTHO/*.ali --config=config-GENERA.yaml --verbosity=5 --threads=$cpu 2> log
        cd ORTHO/
        sed -i -e 's/|/@/g' *GENERA.ali
        find *GENERA.ali > list
        #sed -i -e 's/-BMC_OG-GENERA.ali//g' list
        sed -i -e 's/.ali//g' list
        for f in `cat list`; do ali2fasta.pl \$f.ali; mv \$f.fasta \$f.faa; done
        cd ../
        mkdir ENRICHED
        mv ORTHO/*GENERA.faa ENRICHED/
        echo "GENERA info: 42, 42 done" >> GENERA-SGC.log
        """
    }
    else {
        println "GENERA info: specific enrichment not activated"
        """
        mkdir ENRICHED
        echo "GENERA info: specific enrichment not activated" >> ENRICHED/info
        echo "GENERA info: specific enrichment not activated" >> GENERA-SGC.log
        """
    }
}

//Orthologous enrichment check part
process enrichmentcheck {
	//informations

	//input output
    input:
    file 'ENRICHED' from enriched_ch
    file 'list' from specificlist_ch3
    val confirmog from params.confirmog
    file "GENERA-SGC.log" from log_ch8

    output:
    file 'SPECIFIC' into specenrich_ch
    file 'verified-specific.list' into enrichedlist_ch
    file "GENERA-SGC.log" into log_ch9

    //script
    script:
    if (params.specific == 'yes') {
        println "GENERA info: running specific enrichment ckeck"
        """
        cd ENRICHED
        $confirmog ../list
        cd ../
        mkdir SPECIFIC
        cp ENRICHED/verified-specific.list .
        for f in `cat verified-specific.list`; do mv ENRICHED/\$f SPECIFIC/; done
        """
    }
    else {
        println "GENERA info: running specific enrichment ckeck not activated"
        """
        mkdir SPECIFIC
        echo "GENERA info: specific enrichment check not activated" >> verified-specific.list
        echo "GENERA info: specific enrichment check not activated" >> SPECIFIC/info
        echo "GENERA info: specific enrichment check not activated" >> GENERA-SGC.log
        """
    }
}

//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'core-OG.list' from coreOGlist_ch
    file 'CORE' from core_ch
    file 'SPECIFIC' from specenrich_ch
    file 'verified-specific.list' from enrichedlist_ch
    file 'ANVIO-files' from anviofiles_ch
    file 'OG-BMC' from ogbmc_ch
    file 'Redo' from redo1_ch
    val version from params.version
    file "GENERA-SGC.log" from log_ch9

    output:
    file 'coreGenes' into coreFINAL_ch
    file 'specificGenes' into specificFINAL_ch
    file 'ANVIO_pangenomic-files' into specificfiltFINAL_ch
    file 'final-specific.list' into finalspecific_ch
    file 'Orthologous' into orthoSeqFINAL_ch2
    file 'REDO' into redoFINAL_ch
    file "GENERA-Orthology.log" into logFINAL_ch

    //script
    script:
    """
    #core
    mkdir coreGenes
    mv CORE/* coreGenes/
    #specific
    mkdir specificGenes
    mv SPECIFIC/* specificGenes
    mv ANVIO-files ANVIO_pangenomic-files
    mv verified-specific.list final-specific.list
    #Orthologous
    mkdir Orthologous
    mv OG-BMC/* Orthologous/
    #Redo
    mkdir REDO
    mv Redo/* REDO/
    #Log
    cp GENERA-SGC.log GENERA-Orthology.log
    echo VERSION: >> GENERA-Orthology.log
    echo $version >> GENERA-Orthology.log
    """
}
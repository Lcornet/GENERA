#!/usr/bin/env nextflow

/*
========================================================================================
                        OGsEnrichment
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

    Version: 1.0.1 

    Usage:
    
    The typical command for running the pipeline is as follows:

    OGsEnrichment.nf ... 
    
    Mandatory arguments:
    --OG                     Path to OG directory in fasta format 
    --representative         Path to representative proteomes in fasta format
    --org                    Path to organism to add in fasta format, prot of DNA

    Optional arguments:
    --orgType                Specify org type, DNA or prot, default = prot
    --align                  Specify alignment of OGs, yes or no, default = yes
    --ftaligner              Specify forty-two aligner mode, blast or exonerate, default = blast
    --degap                  Degap OGs affter enrichment and realign (only if --org is prot), yes or no, default = no
    --verbosity              Specify the verbosity of leel, default = 5
    --taxdump 		         Path to taxdump directory - automatic setup by default
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
params.OG = null
if (params.OG == null) {
	exit 1, "Path to OG directory in fasta format"
}

//Path to representative : Mandatory
params.representative = null
if (params.representative == null) {
	exit 1, "Path to representative genomes/proteomes in fasta format"
}

//Path to organisms : Mandatory
params.org = null
if (params.org == null) {
	exit 1, "Path to organism to add genomes/proteomes in fasta format"
}

//OrGType                
params.orgType = 'prot'

//mode
params.verbosity = '5'

//align 
params.align = 'yes'

//42aligner
params.ftaligner = 'blast'

//degap
params.degap = 'no'

//cpu
params.cpu = '1'

//outdir
params.outdir='GENERA-OGsEnrichment'

//Path to taxdump
params.taxdump = 'local'

//Not specified by user
//Path to project dir taxdump
taxdir = "$workflow.projectDir" + '/taxdump'
workingdir = file(taxdir)

//Path to companion
params.companion = '/opt/companion/OGsEnrichment_companion.py'

//version
params.version = '1.0.1'

/*
CORE PROGRAM
*/

//Load input files
taxdump_ch = Channel.fromPath(params.taxdump)
og_ch = Channel.fromPath(params.OG)
representative_ch = Channel.fromPath(params.representative)
org_ch = Channel.fromPath(params.org)


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

//OGs alignment
process alignment {
	//informations

	//input output
    input:
    file '*.f*' from og_ch

    output:
    file 'aligned/*' into ali_ch1
    file "GENERA-OGsEnrichment.log" into log_ch1

    //script
    script:
    println "GENERA info: OGs alignment"
    if (params.align != 'yes'){
        println "GENERA info: OGs alignment"
        """
        mkdir aligned
        mkdir OGs
        cp .f/*.faa OGs
        cd OGs/
        find *.faa | cut -f1 -d"." > list
        mv list ../
        cd ../
        for f in `cat list`; do muscle3.8.31_i86linux64 -in OGs/\$f.faa -out aligned/\$f.faa; done    
        echo "GENERA info: OGs alignment" >> GENERA-OGsEnrichment.log
        """
    }
    else {
        println "GENERA info: alignment not activated"
        """
        mkdir aligned
        cp .f/*.faa aligned/
        echo "GENERA info: alignment not activated" >> GENERA-OGsEnrichment.log
        """
    }
}

//DNA translation
process Enrichment {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file '*' from org_ch
    file '*' from representative_ch
    file 'aligned/*' from ali_ch1
    val taxdir from taxdir_ch1
    val companion from params.companion
    val verbosity from params.verbosity
    val aligner from params.ftaligner
    val cpu from params.cpu
    file "GENERA-OGsEnrichment.log" from log_ch1

    output:
    file 'ENRICHED/*' into enriched_ch
    file "GENERA-OGsEnrichment.log" into log_final

    //script
    script:
    if (params.orgType == 'prot'){
        println "GENERA info: org type is set as prot" 
        if (params.degap == 'no'){
            println "GENERA info: degap not activated"
            """
            echo "GENERA info: degap not activated" >> GENERA-OGsEnrichment.log
            echo "GENERA info: 42, org type is set as prot" >> GENERA-OGsEnrichment.log
            #Reference organism part
            mkdir ref-banks
            cp representative/*.faa ref-banks
            cd ref-banks/
            find *.faa > list
            sed -i -e 's/.faa//g' list
            for f in `cat list`; do inst-abbr-ids.pl \$f.faa --id-regex=:DEF; mv \$f-abbr.faa \$f.faa; sed -i -e 's/|//g' \$f.faa; done
            for f in `cat list`; do makeblastdb -in \$f.faa -dbtype prot -parse_seqids -out \$f; done
            #for f in `cat list `; do echo ".psq" ; done > end.list
            #paste list end.list > part1
            #sed -i -e 's/\t//g' part1
            #paste list part1 > ref-bank-mapper.idm
            paste list list > ref-bank-mapper.idm
            cd ../
            echo "GENERA info: 42, representative bank done" >> GENERA-OGsEnrichment.log

            #Define Queries
            cp ref-banks/list .
            $companion list --mode=queries
            echo "GENERA info: 42, queries done" >> GENERA-OGsEnrichment.log

            #Part for org to add
            mkdir org-to-add
            cp org/*.faa org-to-add/
            cd org-to-add/
            find *.faa > list
            sed -i -e 's/.faa//g' list
            for f in `cat list`; do makeblastdb -in \$f.faa -dbtype prot -parse_seqids -out \$f; done
            #for f in `cat list`; do echo ".psq"; done  > end.list
            #paste list end.list > part1
            #sed -i -e 's/\t//g' part1
            #paste list part1 > bank-mapper.idm
            paste list list > bank-mapper.idm
            cd ..
            echo "GENERA info: 42, org bank done" >> GENERA-OGsEnrichment.log

            #Part OGs
            mkdir ORTHO
            cp aligned/*.faa ORTHO/
            cd ORTHO/
            find *.faa > list
            sed -i -e 's/.faa//g' list
            for f in `cat list`; do fasta2ali.pl \$f.faa; $companion \$f.ali --mode=format; done
            cd ../
            echo "GENERA info: 42, OGs done" >> GENERA-OGsEnrichment.log

            #Generate yaml
            yaml-generator-42.pl --run_mode=phylogenomic --out_suffix=-GENERA --queries queries.idl --evalue=1e-05 --homologues_seg=yes \
            --max_target_seqs=10000 --templates_seg=no --bank_dir org-to-add --bank_suffix=.psq --bank_mapper org-to-add/bank-mapper.idm --ref_brh=on \
            --ref_bank_dir ref-banks --ref_bank_suffix=.psq --ref_bank_mapper ref-banks/ref-bank-mapper.idm --ref_org_mul=0.3 --ref_score_mul=0.99 \
            --trim_homologues=on --ali_keep_lengthened_seqs=keep --aligner_mode=$aligner --tax_reports=off --tax_dir $taxdir \
            --megan_like --tol_check=off
            echo "GENERA info: 42, Yaml done" >> GENERA-OGsEnrichment.log

            #run forty-two
            forty-two.pl ORTHO/*-sp.ali --config=config-GENERA.yaml --verbosity=$verbosity --threads=$cpu 2> log
            cd ORTHO/
            sed -i -e 's/|/@/g' *GENERA.ali
            find *GENERA.ali > list
            sed -i -e 's/-sp-GENERA.ali//g' list
            for f in `cat list`; do $companion \$f-sp-GENERA.ali --mode=restore; ali2fasta.pl \$f-GENERA.ali; mv \$f-GENERA.fasta \$f-GENERA.faa; done
            cd ../
            mkdir ENRICHED
            mv ORTHO/*GENERA.faa ENRICHED/
            echo "GENERA info: 42, 42 done" >> GENERA-OGsEnrichment.log
            """
        }
        else {
            println "GENERA info: degap activated"
            """
            echo "GENERA info: degap activated" >> GENERA-OGsEnrichment.log
            echo "GENERA info: 42, org type is set as prot" >> GENERA-OGsEnrichment.log
            #Reference organism part
            mkdir ref-banks
            cp representative/*.faa ref-banks
            cd ref-banks/
            find *.faa > list
            sed -i -e 's/.faa//g' list
            for f in `cat list`; do inst-abbr-ids.pl \$f.faa --id-regex=:DEF; mv \$f-abbr.faa \$f.faa; sed -i -e 's/|//g' \$f.faa; done
            for f in `cat list`; do makeblastdb -in \$f.faa -dbtype prot -parse_seqids -out \$f; done
            #for f in `cat list `; do echo ".psq" ; done > end.list
            #paste list end.list > part1
            #sed -i -e 's/\t//g' part1
            #paste list part1 > ref-bank-mapper.idm
            paste list list > ref-bank-mapper.idm
            cd ../
            echo "GENERA info: 42, representative bank done" >> GENERA-OGsEnrichment.log

            #Define Queries
            cp ref-banks/list .
            $companion list --mode=queries
            echo "GENERA info: 42, queries done" >> GENERA-OGsEnrichment.log

            #Part for org to add
            mkdir org-to-add
            cp org/*.faa org-to-add/
            cd org-to-add/
            find *.faa > list
            sed -i -e 's/.faa//g' list
            for f in `cat list`; do makeblastdb -in \$f.faa -dbtype prot -parse_seqids -out \$f; done
            #for f in `cat list`; do echo ".psq"; done  > end.list
            #paste list end.list > part1
            #sed -i -e 's/\t//g' part1
            #paste list part1 > bank-mapper.idm
            paste list list > bank-mapper.idm
            cd ..
            echo "GENERA info: 42, org bank done" >> GENERA-OGsEnrichment.log

            #Part OGs
            mkdir ORTHO
            cp aligned/*.faa ORTHO/
            cd ORTHO/
            find *.faa > list
            sed -i -e 's/.faa//g' list
            for f in `cat list`; do fasta2ali.pl \$f.faa; $companion \$f.ali --mode=format; done
            cd ../
            echo "GENERA info: 42, OGs done" >> GENERA-OGsEnrichment.log

            #Generate yaml
            yaml-generator-42.pl --run_mode=phylogenomic --out_suffix=-GENERA --queries queries.idl --evalue=1e-05 --homologues_seg=yes \
            --max_target_seqs=10000 --templates_seg=no --bank_dir org-to-add --bank_suffix=.psq --bank_mapper org-to-add/bank-mapper.idm --ref_brh=on \
            --ref_bank_dir ref-banks --ref_bank_suffix=.psq --ref_bank_mapper ref-banks/ref-bank-mapper.idm --ref_org_mul=0.3 --ref_score_mul=0.99 \
            --trim_homologues=on --ali_keep_lengthened_seqs=keep --aligner_mode=$aligner --tax_reports=off --tax_dir $taxdir \
            --megan_like --tol_check=off
            echo "GENERA info: 42, Yaml done" >> GENERA-OGsEnrichment.log

            #run forty-two
            forty-two.pl ORTHO/*-sp.ali --config=config-GENERA.yaml --verbosity=$verbosity --threads=$cpu 2> log
            cd ORTHO/
            sed -i -e 's/|/@/g' *GENERA.ali
            find *GENERA.ali > list
            sed -i -e 's/-sp-GENERA.ali//g' list
            for f in `cat list`; do $companion \$f-sp-GENERA.ali --mode=restore; ali2fasta.pl \$f-GENERA.ali --degap; mv \$f-GENERA.fasta \$f-GENERA.faa; done
            #realign
            mkdir REALIGNED
            mkdir OGs
            mv *GENERA.faa OGs
            cd OGs/
            find *.faa | cut -f1 -d"." > list
            mv list ../
            cd ../
            for f in `cat list`; do muscle3.8.31_i86linux64 -in OGs/\$f.faa -out REALIGNED/\$f.faa; done
            cd ../
            mkdir ENRICHED
            mv ORTHO/REALIGNED/*.faa ENRICHED/
            echo "GENERA info: 42, 42 done" >> GENERA-OGsEnrichment.log
            """
        }
    }
    else if (params.orgType == 'DNA') {
        println "GENERA info: rog type is set as DNA" 
        """
        echo "GENERA info: 42, org type is set as DNA" >> GENERA-OGsEnrichment.log
        #Reference organism part
        mkdir ref-banks
        cp representative/*.faa ref-banks
        cd ref-banks/
        find *.faa > list
        sed -i -e 's/.faa//g' list
        for f in `cat list`; do inst-abbr-ids.pl \$f.faa --id-regex=:DEF; mv \$f-abbr.faa \$f.faa; sed -i -e 's/|//g' \$f.faa; done
        for f in `cat list`; do makeblastdb -in \$f.faa -dbtype prot -parse_seqids -out \$f; done
        #for f in `cat list `; do echo ".psq" ; done > end.list
        #paste list end.list > part1
        #sed -i -e 's/\t//g' part1
        #paste list part1 > ref-bank-mapper.idm
        paste list list > ref-bank-mapper.idm
        cd ../
        echo "GENERA info: 42, representative bank done" >> GENERA-OGsEnrichment.log

        #Define Queries
        cp ref-banks/list .
        $companion list --mode=queries
        echo "GENERA info: 42, queries done" >> GENERA-OGsEnrichment.log

        #Part for org to add
        mkdir org-to-add
        cp org/*.fna org-to-add/
        cd org-to-add/
        find *.fna > list
        sed -i -e 's/.fna//g' list
        for f in `cat list`; do makeblastdb -in \$f.fna -dbtype nucl -parse_seqids -out \$f; done
        #for f in `cat list`; do echo ".nsq"; done  > end.list
        #paste list end.list > part1
        #sed -i -e 's/\t//g' part1
        #paste list part1 > bank-mapper.idm
        paste list list > bank-mapper.idm
        cd ..
        echo "GENERA info: 42, org bank done" >> GENERA-OGsEnrichment.log

        #Part OGs
        mkdir ORTHO
        cp aligned/*.faa ORTHO/
        cd ORTHO/
        find *.faa > list
        sed -i -e 's/.faa//g' list
        for f in `cat list`; do fasta2ali.pl \$f.faa; $companion \$f.ali --mode=format; done
        cd ../
        echo "GENERA info: 42, OGs done" >> GENERA-OGsEnrichment.log

        #Generate yaml
        yaml-generator-42.pl --run_mode=phylogenomic --out_suffix=-GENERA --queries queries.idl --evalue=1e-05 --homologues_seg=yes \
        --max_target_seqs=10000 --templates_seg=no --bank_dir org-to-add --bank_suffix=.nsq --bank_mapper org-to-add/bank-mapper.idm --ref_brh=on \
        --ref_bank_dir ref-banks --ref_bank_suffix=.psq --ref_bank_mapper ref-banks/ref-bank-mapper.idm --ref_org_mul=0.3 --ref_score_mul=0.99 \
        --trim_homologues=on --ali_keep_lengthened_seqs=keep --aligner_mode=$aligner --tax_reports=off --tax_dir $taxdir \
        --megan_like --tol_check=off
        echo "GENERA info: 42, Yaml done" >> GENERA-OGsEnrichment.log

        #run forty-two
        forty-two.pl ORTHO/*-sp.ali --config=config-GENERA.yaml --verbosity=$verbosity --threads=$cpu 2> log
        cd ORTHO/
        sed -i -e 's/|/@/g' *GENERA.ali
        find *GENERA.ali > list
        sed -i -e 's/-sp-GENERA.ali//g' list
        for f in `cat list`; do $companion \$f-sp-GENERA.ali --mode=restore; ali2fasta.pl \$f-GENERA.ali; mv \$f-GENERA.fasta \$f-GENERA.faa; done
        cd ../
        mkdir ENRICHED
        mv ORTHO/*GENERA.faa ENRICHED/
        echo "GENERA info: 42, 42 done" >> GENERA-OGsEnrichment.log 
        """
    }
}
#!/usr/bin/env nextflow

/*
========================================================================================
                        GENERA contams
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

    nextflow .... 
    
    Mandatory arguments:
    --genomes                Specify the directory with genomes, shoudl be 'GENERA-input/'

    Optional arguments:
    --mode                   Specify mode, all or specific tool, default = all
    --taxlevel               Specify taxonomic level to use for Physeter/Kraken, default = phylum
    --automode               Specify methods for main organism detection in Physeter/Kraken, label_first or count_first, default = label_first
    --idl                    Path to idl file for Physeter, automatic build by default according to taxlevel option
    --ckcompleteness         Final list: Minimum CheckM completeness, default = 95
    --ckcontamination        Final list: Maximum CheckM contamination, default = 5
    --gunccss                Final list: Maximum GUNC css, default = 0.01
    --guncrrs                Final list: Minimum GUNC rrs, default = 0.5
    --physetercontamination  Final list: Maximum Physeter contamination, default = 100 (unactivated by default)
    --krakencontamination    Final list: Maximum Kraken contamination, default = 100 (unactivated by default)
    --bucompleteness         Final list: Minimum BUSCO completeness, default = 0 (unactivated by default)
    --budups                 Final list: Maximum BUSCO duplication, default = 100 (unactivated by default)
    --numcontigs             Final list: Maximum Number of contigs, default = 1000
    --ext                    Specify the extention of genomes, fna or fa or fasta, default = fna
    --dbdir                  Path to GENERA-DB directory, automatic setup by default
    --taxdump                Path to taxdump, automatic setup by default
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
params.genomes = null
if (params.genomes == null) {
	exit 1, "Specify directory with genomes"
}

//mode
params.mode = 'all'

//ext
params.ext = 'fna'

//DB dir
params.dbdir = 'local'
//Path to db dir
dbdirectory = "$workflow.projectDir" + '/GENERA-DB-bak'
workingdb = file(dbdirectory)

//Path to taxdump
params.taxdump = 'local'

//Not specified by user
//Path to project dir taxdump
taxdir = "$workflow.projectDir" + '/taxdump'
workingdir = file(taxdir)

//automode
params.automode = 'label_first'

//ckcompleteness
params.ckcompleteness = '95'

//ckcontamination
params.ckcontamination = '5'

//gunccss
params.gunccss = '0.01'

//guncrrs
params.guncrrs = '0.5'

//physetercontamination
params.physetercontamination = '100'

//krakencontamination
params.krakencontamination = '100'

//bucompleteness
params.bucompleteness = '0'

//budups
params.budups = '100'

//numcontigs
params.numcontigs = '1000'

//cpu
params.cpu = '1'

//taxo level
params.taxlevel = 'phylum'

//IDL
params.idl = '/opt/companion/contam-labels.idl'

//outdir
params.outdir='GENERA-contams'

//Path to companion
params.companion = '/opt/companion/Contams_companion.py'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
taxdump_ch = Channel.fromPath(params.taxdump)
genomes_ch = Channel.fromPath(params.genomes)
dbdir_ch1 = Channel.fromPath(params.dbdir)
dbdir_ch2 = Channel.fromPath(params.dbdir)
dbdir_ch3 = Channel.fromPath(params.dbdir)


//NCBI Taxonomy, set taxdump if not specifed
process taxonomy {
	//informations

	//input output
    input:
    val taxdump from taxdump_ch
    
    output:
    file "taxdump_path.txt" into taxdump_path1
    val taxdir into taxdir_ch1
    val taxdir into taxdir_ch2

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


//Download DB
//Get gunc bd files
process DBSetUp {
	//informations

	//input output
    input:
    val dbdir from dbdir_ch1
    
    output:
    val databases into dabatses_ch1
    val databases into dabatses_ch2
    val databases into dabatses_ch3
    val databases into dabatses_ch4
    file "db_path.txt" into db_path1
    file 'GENERA-contams.log' into log_ch1

    //script
    script:

    databases = 'na'

    if (params.dbdir == 'local'){
        println "GENERA-INFO: DB directory not specified -> project dir"

        if( !workingdb.exists() ) {
            println "GENERA-INFO: DB directory not found in project dir -> Created"
            if( !workingdb.mkdirs() )    {
                exit 1, "Cannot create working directory"
            }

            databases = workingdb

            """
            echo "GENERA info: DB directory not specified" >> GENERA-contams.log
            echo "GENERA info: DB directory not found in project dir -> Created" >> GENERA-contams.log
            #Set up GUNC DB
            mkdir $workingdb/guncDB
            gunc download_db $workingdb/guncDB

            #Set up Physeter DB
            wget https://figshare.com/ndownloader/files/32405687 -O Cornet-Baurain.tgz
            tar -xzf Cornet-Baurain.tgz
            mkdir $workingdb/Physeter
            mv Cornet-2022-GBIO-Figshare/life-tqmd-of73* $workingdb/Physeter/
            mv Cornet-2022-GBIO-Figshare/contam-labels.idl $workingdb/Physeter/

            #Set up BUSCO DB
            busco --download all
            mv busco_downloads/ $workingdb/

            #kraken
            #kraken2-build --download-library bacteria --db Kbact
            #kraken2-build --build --db Kbact --threads=20

            #DB done
            echo "GENERA info: DB done" >> GENERA-contams.log
            echo $workingdb > db_path.txt
            """
        }
        else {
            println "GENERA-INFO: DB directory found in project dir -> Used"

            databases = workingdb 

 	        """
            echo $workingdb > db_path.txt
            echo "GENERA info: DB directory not specified" >> GENERA-contams.log
            echo "GENERA-INFO: DB directory found in project dir -> Used" >> GENERA-contams.log
		    """           

        }
    }

	else{
        println "GENERA-INFO: DB directory specifies"

        databases = dbdir

		"""
        echo $dbdir > db_path.txt
        echo "GENERA-INFO: DB directory specified" >> GENERA-contams.log
		"""		
    }
}


//format
process format {
	//informations

	//input output
    input:
    file '*' from genomes_ch
    file 'GENERA-contams.log' from log_ch1

    output:
    file 'GENOMES/' into genomes_ch1
    file 'GENOMES/' into genomes_ch2
    file 'GENOMES/' into genomes_ch3
    file 'GENOMES/' into genomes_ch4
    file 'GENOMES/' into genomes_ch5
    file 'SPLIT/' into split_ch1
    file 'SPLIT/' into split_ch2
    file 'GENERA-contams.log' into log_ch2

    //script
    script:
    if (params.ext == 'fna') {
        """
        mkdir TEMP
        cp GENERA-input/*.fna TEMP
        cd TEMP/
        find *.fna > list
        sed -i -e 's/.fna//g' list
        for f in `cat list`; do inst-abbr-ids.pl \$f.fna --id-regex=:DEF --id-prefix=\$f; done
        for f in `cat list`; do inst-split-seqs.pl \$f-abbr.fna --out=-split; done
        cd ../
        mkdir SPLIT
        mv TEMP/*split* SPLIT/
        mkdir GENOMES/
        mv TEMP/*abbr* GENOMES/
        echo "GENERA info: format .fna files" >> GENERA-contams.log
        """
    }
    else if (params.ext == 'fasta') {
        """
        mkdir TEMP
        cp GENERA-input/*.fasta TEMP
        cd TEMP
        find *.fasta > list
        sed -i -e 's/.fasta//g' list
        for f in `cat list`; do mv \$f.fasta \$f.fna; done
        for f in `cat list`; do inst-abbr-ids.pl \$f.fna --id-regex=:DEF --id-prefix=\$f; done
        for f in `cat list`; do inst-split-seqs.pl \$f-abbr.fna --out=-split; done
        cd ../
        mkdir SPLIT
        mv TEMP/*split* SPLIT/
        mkdir GENOMES/
        mv TEMP/*abbr* GENOMES/
        echo "GENERA info: format .fasta files" >> GENERA-contams.log
        """
    }
    else if (params.ext == 'fa') {
        """
        mkdir TEMP
        cp GENERA-input/*.fa TEMP
        cd TEMP
        find *.fa > list
        sed -i -e 's/.fa//g' list
        for f in `cat list`; do mv \$f.fa \$f.fna; done
        for f in `cat list`; do inst-abbr-ids.pl \$f.fna --id-regex=:DEF --id-prefix=\$f; done
        for f in `cat list`; do inst-split-seqs.pl \$f-abbr.fna --out=-split; done
        cd ../
        mkdir SPLIT
        mv TEMP/*split* SPLIT/
        mkdir GENOMES/
        mv TEMP/*abbr* GENOMES/
        echo "GENERA info: format .fa files" >> GENERA-contams.log
        """
    }
}

//Checkm
process checkm {
	//informations

	//input output
    input:
    file '*' from genomes_ch1
    val cpu from params.cpu
    file 'GENERA-contams.log' from log_ch2

    output:
    file 'results.Checkm' into checkm_ch1
    file 'GENERA-contams.log' into log_ch3

    //script
    script:
    if (params.mode == 'checkm') {
        """
        #mkdir GEN
        #cp GENOMES/*abbr.fna GEN/
        checkm lineage_wf -t $cpu -x fna GENOMES runc > checkm.result
        echo "#genome,completeness,contamination,str-hetero" > part1
        tr -s " " < checkm.result | grep "abbr" | cut -f2,14,15,16 -d" " > part2
        sed -i -e 's/ /,/g' part2
        cat part1 part2 > results.Checkm
        echo "GENERA info: run checkm" >> GENERA-contams.log
        """
    }
    else if (params.mode == 'all') {
        """
        #mkdir GEN
        #cp GENOMES/*abbr.fna GEN/
        checkm lineage_wf -t $cpu -x fna GENOMES runc > checkm.result
        echo "#genome,completeness,contamination,str-hetero" > part1
        tr -s " " < checkm.result | grep "abbr" | cut -f2,14,15,16 -d" " > part2
        sed -i -e 's/ /,/g' part2
        cat part1 part2 > results.Checkm
        echo "GENERA info: run checkm" >> GENERA-contams.log
        """
    }
    else {
        """
        echo 'FALSE' > results.Checkm
        echo "GENERA info: checkm not activated" >> GENERA-contams.log
        """
    }
}

//GUNC
process gunc {
	//informations

	//input output
    input:
    file '*' from genomes_ch2
    file 'GENERA-contams.log' from log_ch3
    val databases from dabatses_ch1
    val cpu from params.cpu

    output:
    file 'GUNC/GUNC.progenomes_2.1.maxCSS_level.tsv' into gunc_ch1
    file 'GENERA-contams.log' into log_ch4

    //script
    script:
    if (params.mode == 'gunc') {
        """
        #mkdir GEN
        #cp GENOMES/*abbr.fna GEN/
        mkdir GUNC
        gunc run --db $databases/guncDB/gunc_db_progenomes2.1.dmnd --input_dir GENOMES/ --file_suffix .fna --threads 20 --out_dir GUNC/
        echo "GENERA info: run gunc" >> GENERA-contams.log
        """
    }
    else if (params.mode == 'all') {
        """
        #mkdir GEN
        #cp GENOMES/*abbr.fna GEN/
        mkdir GUNC
        gunc run --db $databases/guncDB/gunc_db_progenomes2.1.dmnd --input_dir GENOMES/ --file_suffix .fna --threads 20 --out_dir GUNC/
        echo "GENERA info: run gunc" >> GENERA-contams.log
        """
    }
    else {
        """
        mkdir GUNC
        cd GUNC
        echo FALSE > GUNC.progenomes_2.1.maxCSS_level.tsv
        cd ../
        echo "GENERA info: gunc not activated" >> GENERA-contams.log
        """
    }
}

//BUSCO
process busco {
	//informations

	//input output
    input:
    file '*' from genomes_ch3
    file 'GENERA-contams.log' from log_ch4
    val databases from dabatses_ch2
    val cpu from params.cpu

    output:
    file 'BUSCO/batch_summary.txt' into busco_ch1
    file 'GENERA-contams.log' into log_ch5

    //script
    script:
    if (params.mode == 'busco') {
        """
        mkdir GEN
        cp GENOMES/*abbr.fna GEN/
        cd GEN/
        find *.fna > list
        sed -i -e 's/.fna//g' list
        for f in `cat list`; do  sed -i -e 's/|/@/g' \$f.fna ; done
        cd ../
        busco -m genome -i GEN/ -o BUSCO --auto-lineage --download_path $databases/busco_downloads/ --cpu $cpu
        echo "GENERA info: run busco" >> GENERA-contams.log
        """
    }
    else if (params.mode == 'all') {
        """
        #mkdir GEN
        #cp GENOMES/*abbr.fna GEN/
        #cd GEN/
        #find *.fna > list
        #sed -i -e 's/.fna//g' list
        #for f in `cat list`; do  sed -i -e 's/|/@/g' \$f.fna ; done
        #cd ../
        #busco -m genome -i GEN/ -o BUSCO --auto-lineage --download_path $databases/busco_downloads/ --cpu $cpu
        #echo "GENERA info: run busco" >> GENERA-contams.log
        mkdir BUSCO/
        cd BUSCO/
        echo 'FALSE' > batch_summary.txt
        cd ../
        echo "GENERA info: busco not activated" >> GENERA-contams.log
        """
    }
    else {
        """
        mkdir BUSCO/
        cd BUSCO/
        echo 'FALSE' > batch_summary.txt
        cd ../
        echo "GENERA info: busco not activated" >> GENERA-contams.log
        """
    }
}

//Physeter
process physeter {
	//informations

	//input output
    input:
    file '*' from split_ch1
    val databases from dabatses_ch3
    val cpu from params.cpu
    val idl from params.idl
    val taxdir from taxdir_ch1
    val taxlevel from params.taxlevel
    val automode from params.automode
    file 'GENERA-contams.log' from log_ch5

    output:
    file 'Physeter.report' into physeter_ch1
    file 'GENERA-contams.log' into log_ch6

    //script
    script:
    if (params.mode == 'physeter') {
        """
        #labeller
        grep genus $taxdir/nodes.dmp | cut -f1 > genus.taxid
        /opt/create-labeler.pl genus.taxid --taxdir=$taxdir --level=$taxlevel --kingdoms=Bacteria Archaea Eukaryota > file.idl

        mkdir GEN
        cp SPLIT/*split.fna GEN/
        cd GEN/
        mv ../file.idl .
        find *.fna > list
        sed -i -e 's/.fna//g' list
        for f in `cat list`; do mv \$f.fna \$f.fasta; done
        #Run Diamond
        for f in `cat list`; do mkdir temp; diamond blastx -d $databases/Physeter/life-tqmd-of73.dmnd \
        -q \$f.fasta -o \$f.blastx -t temp -k 50 -e 1e-10 -f tab -p $cpu ; rm -rf temp; done
        #Run Physeter
        for f in `cat list`; do physeter.pl \$f.blastx --fasta-dir=./ --outfile=\$f.report --taxdir=$taxdir \
        --taxon-list=file.idl --auto-detect --kraken; done
        #--taxon-list=$idl --auto-detect --kraken; done
        #Run kraken-parser.pl
        for f in `cat list`; do /opt/kraken-parser.pl \$f-kraken.tsv --taxdir=$taxdir --outfile=\$f-parsed.report \
        --taxon-list=file.idl --auto-detect=$automode; done
        cat *parsed.report > ../Physeter.report
        cd ../
        echo "GENERA info: run Physeter" >> GENERA-contams.log
        """
    }
    else if (params.mode == 'all') {
        """
        #labeller
        grep genus $taxdir/nodes.dmp | cut -f1 > genus.taxid
        /opt/create-labeler.pl genus.taxid --taxdir=$taxdir --level=$taxlevel --kingdoms=Bacteria Archaea Eukaryota > file.idl

        mkdir GEN
        cp SPLIT/*split.fna GEN/
        cd GEN/
        mv ../file.idl .
        find *.fna > list
        sed -i -e 's/.fna//g' list
        for f in `cat list`; do mv \$f.fna \$f.fasta; done
        #Run Diamond
        for f in `cat list`; do mkdir temp; diamond blastx -d $databases/Physeter/life-tqmd-of73.dmnd \
        -q \$f.fasta -o \$f.blastx -t temp -k 50 -e 1e-10 -f tab -p $cpu ; rm -rf temp; done
        #Run Physeter
        for f in `cat list`; do physeter.pl \$f.blastx --fasta-dir=./ --outfile=\$f.report --taxdir=$taxdir \
        --taxon-list=file.idl --auto-detect --kraken; done
        #--taxon-list=$idl --auto-detect --kraken; done
        #Run kraken-parser.pl
        for f in `cat list`; do /opt/kraken-parser.pl \$f-kraken.tsv --taxdir=$taxdir --outfile=\$f-parsed.report \
        --taxon-list=file.idl --auto-detect=$automode; done
        cat *parsed.report > ../Physeter.report
        cd ../
        echo "GENERA info: run Physeter" >> GENERA-contams.log
        """
    }
    else {
        """
        echo 'FALSE' > Physeter.report
        echo "GENERA info: Physeter not activated" >> GENERA-contams.log
        """
    }
}

//kraken
process kraken {
	//informations

	//input output
    input:
    file '*' from split_ch2
    val databases from dabatses_ch4
    val cpu from params.cpu
    val taxdir from taxdir_ch2
    val taxlevel from params.taxlevel
    val automode from params.automode
    file 'GENERA-contams.log' from log_ch6

    output:
    file 'KRAKEN/' into kraken_ch1
    file 'Kraken.report' into kraken_ch2
    file 'GENERA-contams.log' into log_ch7

    //script
    script:
    if (params.mode == 'kraken') {
        """
        #labeller
        grep genus $taxdir/nodes.dmp | cut -f1 > genus.taxid
        /opt/create-labeler.pl genus.taxid --taxdir=$taxdir --level=$taxlevel --kingdoms=Bacteria Archaea Eukaryota > file.idl

        mkdir GEN
        cp SPLIT/*split.fna GEN/
        cd GEN/
        mv ../file.idl .
        find *split.fna > list
        sed -i -e 's/.fna//g' list
        #Run kraken
        for f in `cat list`; do kraken2 --use-names --db $databases/Kraken/ \$f.fna --threads $cpu \
        --report \$f.report > \$f.kraken; done
        #labeller
        for f in `cat list`; do kraken2 --use-names --db $databases/Kraken/ \$f.fna --threads $cpu \
        --report \$f.report > \$f.kraken; done  
        for f in `cat list`; do /opt/kraken-parser.pl \$f.report --taxdir=$taxdir --outfile=\$f-parsed.report \
        --taxon-list=file.idl --auto-detect=$automode; done
        #store
        cd ../
        mkdir KRAKEN/
        mv GEN/*.report GEN/*.kraken KRAKEN/
        cat KRAKEN/*parsed.report > Kraken.report
        echo "GENERA info: run Kraken2" >> GENERA-contams.log
        """
    }
    else if (params.mode == 'all') {
        """
        #labeller
        grep genus $taxdir/nodes.dmp | cut -f1 > genus.taxid
        /opt/create-labeler.pl genus.taxid --taxdir=$taxdir --level=$taxlevel --kingdoms=Bacteria Archaea Eukaryota > file.idl

        mkdir GEN
        cp SPLIT/*split.fna GEN/
        cd GEN/
        mv ../file.idl .
        find *split.fna > list
        sed -i -e 's/.fna//g' list
        #Run kraken
        for f in `cat list`; do kraken2 --use-names --db $databases/Kraken/ \$f.fna --threads $cpu \
        --report \$f.report > \$f.kraken; done
        #labeller
        for f in `cat list`; do kraken2 --use-names --db $databases/Kraken/ \$f.fna --threads $cpu \
        --report \$f.report > \$f.kraken; done  
        for f in `cat list`; do /opt/kraken-parser.pl \$f.report --taxdir=$taxdir --outfile=\$f-parsed.report \
        --taxon-list=file.idl --auto-detect=$automode; done
        #store
        cd ../
        mkdir KRAKEN/
        mv GEN/*.report GEN/*.kraken KRAKEN/
        cat KRAKEN/*parsed.report > Kraken.report
        echo "GENERA info: run Kraken2" >> GENERA-contams.log
        """
    }
    else {
        """
        echo 'FALSE' > FALSE.kraken
        echo 'FALSE' > FALSE.report
        #store
        mkdir KRAKEN/
        mv *.report *.kraken KRAKEN/
        echo "GENERA info: Kraken2 not activated" > Kraken.report
        echo "GENERA info: Kraken2 not activated" >> GENERA-contams.log
        """
    }
}

//quast
process quast {
	//informations

	//input output
    input:
    file '*' from genomes_ch4
    val cpu from params.cpu
    file 'GENERA-contams.log' from log_ch7

    output:
    file 'quast.report' into quast_ch1
    file 'GENERA-contams.log' into log_ch8

    //script
    script:
    if (params.mode == 'quast') {
        """
        find GENOMES/* > list
        sed -i -e 's/.fna//g' list
        #Run quast
        for f in `cat list`; do quast.py \$f.fna -t $cpu -o QUAST-\$f; done
        cat QUAST-GENOMES/*/transposed_report.tsv > quast.temp
        sed '/^\$/d' quast.temp > temp
        mv temp quast.report
        echo "GENERA info: run quast" >> GENERA-contams.log
        """
    }
    else if (params.mode == 'all') {
        """
        find GENOMES/* > list
        sed -i -e 's/.fna//g' list
        #Run quast
        for f in `cat list`; do quast.py \$f.fna -t $cpu -o QUAST-\$f; done
        cat QUAST-GENOMES/*/transposed_report.tsv > quast.temp
        sed '/^\$/d' quast.temp > temp
        mv temp quast.report
        echo "GENERA info: run quast" >> GENERA-contams.log
        """
    }
    else {
        """
        echo 'FALSE' > quast.report
        echo "GENERA info: quast not activated" >> GENERA-contams.log
        """
    }
}

//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'GENOMES/' from genomes_ch5
    file 'results.Checkm' from checkm_ch1
    file 'GUNC.progenomes_2.1.maxCSS_level.tsv' from gunc_ch1
    file 'batch_summary.txt*' from busco_ch1
    file 'Physeter.report' from physeter_ch1
    file 'KRAKEN/' from kraken_ch1
    file 'Kraken.report' from kraken_ch2
    file 'quast.report' from quast_ch1
    file 'GENERA-contams.log' from log_ch8
    val companion from params.companion
    val ckcompleteness from params.ckcompleteness
    val ckcontamination from params.ckcontamination
    val gunccss from params.gunccss
    val guncrrs from params.guncrrs
    val physetercontamination from params.physetercontamination
    val krakencontamination from params.krakencontamination 
    val bucompleteness from params.bucompleteness
    val budups from params.budups
    val numcontigs from params.numcontigs
    val version from params.version

    output:
    file 'GENERA-contams.table' into finnalTable_ch
    file 'positive-list.txt' into finalpass_ch
    file 'GENERA-contams.log' into finallog_ch
    file 'Checkm.results' into checkmfinal_ch
    file 'GUNC.results' into guncfinal_ch
    file 'Busco.results' into buscofinal_ch
    file 'Physeter.results' into physeterfinal_ch
    file 'KRAKEN-results' into  krakenfinal_ch
    file 'Kraken.results' into  krakenRfinal_ch
    file 'Quast.results' into quastfinal_ch

    //script
    script:
    if (params.mode == 'all') {
        """
	echo "GENERA info: not All mode, no final table or positive list" > GENERA-contams.table
        echo "GENERA info: not All mode, no final table or positive list" > positive-list.txt
        cp results.Checkm Checkm.results
        cp GUNC.progenomes_2.1.maxCSS_level.tsv GUNC.results
        cp batch_summary.txt* Busco.results
        cp Physeter.report Physeter.results
        cp quast.report Quast.results
        cp Kraken.report Kraken.results
        mkdir KRAKEN-results
        cp -r KRAKEN/* KRAKEN-results/
        echo "GENERA info: making final table" >> GENERA-contams.log
        echo "Version:" >> GENERA-contams.log
        echo $version >> GENERA-contams.log
        """
    }
    else {
        """
        echo "GENERA info: not All mode, no final table or positive list" > GENERA-contams.table
        echo "GENERA info: not All mode, no final table or positive list" > positive-list.txt
        cp results.Checkm Checkm.results
        cp GUNC.progenomes_2.1.maxCSS_level.tsv GUNC.results
        cp batch_summary.txt* Busco.results
        cp Physeter.report Physeter.results
        cp quast.report Quast.results
        cp Kraken.report Kraken.results
        mkdir KRAKEN-results
        cp -r KRAKEN/* KRAKEN-results/
        echo "GENERA info: not All mode, no final table or positive list" >> GENERA-contams.log
        echo "Version:" >> GENERA-contams.log
        echo $version >> GENERA-contams.log
        """      
    }
}

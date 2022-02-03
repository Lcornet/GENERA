#!/usr/bin/env nextflow

/*
========================================================================================
                        OGsRtranslate
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

    nextflow run  OGsRtranslate.nf --OG=OGs --bank=bank
    
    Mandatory arguments:
    --OG                     Path to protein OG directory in fasta format
    --bank                   Path to Genomes directory in fasta format

    Optional arguments:
    --verbosity              Specify the verbosity of leel, default = 5

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

//Path to OGs : Mandatory
params.bank = null
if (params.bank == null) {
	exit 1, "Path to Genomes directory in fasta format"
}

//mode
params.verbosity = '5'

//cpu
params.cpu = '1'

//outdir
params.outdir='GENERA-OGsRtranslate'

//Path to companion
params.companion = '/opt/companion/OGsRtranslate_companion.py'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
og_ch = Channel.fromPath(params.OG)
bank_ch = Channel.fromPath(params.bank)

//OGs alignment
process alignment {
	//informations

	//input output
    input:
    file '*.f*' from og_ch

    output:
    file 'aligned/*' into ali_ch1
    file "GENERA-OGsRtranslate.log" into log_ch1

    //script
    script:
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
    echo "GENERA info: OGs alignment" >> GENERA-OGsRtranslate.log
    """
}

//DNA translation
process DNAtranslation {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file '*' from bank_ch
    file '*' from ali_ch1
    val companion from params.companion
    val verbosity from params.verbosity
    val version from params.version
    file "GENERA-OGsRtranslate.log" from log_ch1

    output:
    file 'DNA/*.fna' into leelFINAL_ch
    file "GENERA-OGsRtranslate.log" into logFINAL_ch

    //script
    script:
    println "GENERA info: translate"
    """
    #Banks part
    mkdir banks
    cp bank/*.fna banks/
    cd banks/
    find *.fna > list
    sed -i -e 's/.fna//g' list
    for f in `cat list`; do makeblastdb -in \$f.fna -dbtype nucl -parse_seqids -out \$f; done
    cd ../
    echo "GENERA info: making bank" >> GENERA-OGsRtranslate.log

    #Yaml part
    cd banks/
    $companion list --mode=yaml
    mv yaml.part2 ../
    cd ../
    cat /opt/yaml.part1 yaml.part2 > leel.yaml
    echo "GENERA info: making yaml" >> GENERA-OGsRtranslate.log

    #Leel part
    mkdir ali
    mv *.faa ali/
    cd ali/
    fasta2ali.pl *.faa
    for f in *.ali; do $companion \$f --mode=fsp; done
    cd ../
    leel.pl --config=leel.yaml --verbosity=$verbosity ali/*-sp.ali 2> logleel
    cd ali/
    for f in *-sp-nucl.ali; do $companion \$f --mode=restore; mv \$f \$f.bak; done
    ali2fasta.pl *nucl.ali
    cd ../
    echo "GENERA info: running leel" >> GENERA-OGsRtranslate.log 

    #Publication of results
    mkdir DNA/
    mv ali/*nucl.fasta DNA/
    cd DNA/
    find *.fasta > list
    sed -i -e 's/.fasta//g' list
    for f in `cat list`; do mv \$f.fasta \$f.fna; done
    cd ../
    echo VERSION: >> GENERA-OGsRtranslate.log
    echo $version >> GENERA-OGsRtranslate.log
    """
}
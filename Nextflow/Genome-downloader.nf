#!/usr/bin/env nextflow
/*
========================================================================================
                         Genome Downloader
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
    This tool download genomes from the NCBI.
    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow Genome-downloader.nf --reftaxolevel=phylum --refgroup=Cyanobacteria 
    
    Mandatory arguments:
    --group          Group of interest
    --taxolevel      Taxonomic level of group - Choice between five taxa levels: phylum, class, order, family, genus           

    Optional arguments:
    --taxdump 		    Path to taxdump directory - automatic setup by default
    --genbank           Download GenBank metadata, activated by default - yes or no. Needed if refgenbank or outgenbank is activated
    
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

//Specified by user
//Name of the Refgroup : Mandatory
params.group = null
if (params.group == null) {
	exit 1, "RefSeq group mandatory"
}

//Choice between four taxa levels for refgroup: phylum, class, order, family: Mandatory
params.taxolevel = null
if (params.taxolevel == null) {
	exit 1, "Taxonomic levels: Choice between four taxa levels: phylum, class, order, family, genus "
}

//Path to companion
params.companion = '/opt/COMPANION/Genome-downloader_companion.py'

//Path to taxdump
params.taxdump = 'local'

//GenBank metadata activated by default
params.genbank = 'yes'


//Not specified by user
//Path to project dir taxdump
taxdir = "$workflow.projectDir" + '/taxdump'
workingdir = file(taxdir)

//outdir
params.outdir='Genome-downloader_output'


/*
CORE PROGRAM
*/

//Load input files
taxdump_ch = Channel.fromPath(params.taxdump)

//Taxonomy, set taxdump if not specifed
process Taxonomy {
	//informations

	//input output
    input:
    val taxdump from taxdump_ch
    
    output:
    file "taxdump_path.txt" into taxdump_path1
    file "taxdump_path.txt" into taxdump_path2
    file "taxdump_path.txt" into taxdump_path3
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

//Download RefSeq metadata, compute taxonomy file and prudce download ftp file
process RefSeq {
	//informations

	//input output
    input:
    file taxdump from taxdump_path1
    val companion from params.companion
    
    output:
    file "ftp.sh" into refseq_ftp1
    file "GCF.tax" into refseq_tax1

    //script
    script:

    """
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O refseq_sum.txt
    $companion refseq_sum.txt --mode=sum
    grep -v "#" refseq_sum-filt.txt | cut -f1 > GCF.list
    fetch-tax.pl GCF.list  --taxdir=\$(<taxdump_path.txt) --item-type=taxid --levels=phylum class order family genus
    grep -v "#" refseq_sum-filt.txt | cut -f20 > ftp.list
    grep -v "#" refseq_sum-filt.txt | cut -f20 | cut -f10 -d"/" > names.list
    for f in `cat ftp.list `; do echo "/"; done > slash.list
    for f in `cat ftp.list `; do echo "_genomic.fna.gz"; done > end1.list
    for f in `cat ftp.list `; do echo "wget "; done > get.list
    for f in `cat ftp.list `; do echo " -O "; done > out.list
    for f in `cat ftp.list `; do echo ".fna.gz"; done > end2.list
    cut -f1,2 -d"_" names.list > id.list
    paste get.list ftp.list slash.list names.list end1.list out.list id.list end2.list > ftp.sh
    sed -i -e 's/\t//g' ftp.sh
    """
}

//Download Genbank metadata, compute taxonomy file and prudce download ftp file: OPTIONAL
process GenBank {

    //informations

	//input output
    input:
    file taxdump from taxdump_path2
    val companion from params.companion
    
    output:
    file "GCA-ftp.sh" into genbank_ftp1
    file "GCA.tax" into genbank_tax1

    //script
    script:
    if (params.genbank == 'yes'){
        println "Add GenBank Genomes activated"
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -O genbank_sum.txt
        $companion genbank_sum.txt --mode=sum
        grep -v "#" genbank_sum-filt.txt | cut -f1 > GCA.list
        fetch-tax.pl GCA.list  --taxdir=\$(<taxdump_path.txt) --item-type=taxid --levels=phylum class order family genus
        grep -v "#" genbank_sum-filt.txt | cut -f20 > ftp.list
        grep -v "#" genbank_sum-filt.txt | cut -f20 | cut -f10 -d"/" > names.list
        for f in `cat ftp.list `; do echo "/"; done > slash.list
        for f in `cat ftp.list `; do echo "_genomic.fna.gz"; done > end1.list
        for f in `cat ftp.list `; do echo "wget "; done > get.list
        for f in `cat ftp.list `; do echo " -O "; done > out.list
        for f in `cat ftp.list `; do echo ".fna.gz"; done > end2.list
        cut -f1,2 -d"_" names.list > id.list
        paste get.list ftp.list slash.list names.list end1.list out.list id.list end2.list > GCA-ftp.sh
        sed -i -e 's/\t//g' GCA-ftp.sh
        """
    }
    else {
        println "Add GenBank Genomes NOT activated"
        """
        echo "Add GenBank Genomes NOT activated" > GCA.tax
        echo "Add GenBank Genomes NOT activated" > GCA-ftp.sh
        """

    }

}

//RefGroup part
//Get Refseq genome for the reference group, abbr files
process GetGenomesRefseq {
	//informations

	//input output
    input:
    val group from params.group
    val taxa from params.taxolevel
    val companion from params.companion
    file "ftp.sh" from refseq_ftp1
    file "GCF.tax" from refseq_tax1
    
    output:
    file '*-abbr.fna' into refgenomes_ch1
    file 'reduce-ftp.sh' into reduceRefFtp_ch
    file 'GCF.refgroup.uniq' into refgroupRefseqGC_ch

    //script
    script:

    """
    #Produce list of GCF IDs with reference group and taxa levels
    $companion GCF.tax --mode=fetch --taxa=$taxa --refgroup=$group
    for f in `cat GCF.refgroup.uniq`; do grep \$f ftp.sh; done > reduce-ftp.sh
    bash reduce-ftp.sh
    gunzip *.gz
    find *.fna | cut -f1,2 -d"." > fna.list
    for f in `cat fna.list`; do inst-abbr-ids.pl \$f*.fna --id-regex=:DEF --id-prefix=\$f; done
    """
}

//Get GenBank genome for the reference group, abbr files: OPTIONAL
process GetGenomesGenbank {

    //informations

	//input output
    input:
    val group from params.group
    val taxa from params.taxolevel
    val companion from params.companion
    file "GCA-ftp.sh" from genbank_ftp1
    file "GCA.tax" from genbank_tax1
    file 'GCF.refgroup.uniq' from refgroupRefseqGC_ch
    
    output:
    file '*-abbr.fna' into refgenomesGB_ch1
    file 'GCA-reduce-ftp.sh' into reduceRefFtpGB_ch

    //script
    script:
    if (params.genbank == 'yes'){
        println "Add GenBank Genomes activated"
        """
        #Produce list of GCA IDs with reference group and taxa levels
        $companion GCA.tax --mode=fetch --taxa=$taxa --refgroup=$group
        for f in `cat GCA.refgroup.uniq`; do grep \$f GCA-ftp.sh; done > GCA-reduce-ftp.sh
        bash GCA-reduce-ftp.sh
        gunzip *.gz
        find *.fna | cut -f1,2 -d"." > fna.list
        for f in `cat fna.list`; do inst-abbr-ids.pl \$f*.fna --id-regex=:DEF --id-prefix=\$f; done
        #for fix and proceed , false genbank files
        echo "Add GenBank Genomes activated" > FALSE-abbr.fna
        echo "Add GenBank Genomes activated" > FALSE-GCA-reduce-ftp.sh
        """
    }
    else {
        println "Add GenBank Genomes NOT activated"
        """
        echo "Add GenBank Genomes NOT activated" > FALSE-abbr.fna
        echo "Add GenBank Genomes NOT activated" > GCA-reduce-ftp.sh
        """

    }

}

process CombineGenomes {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file taxdump from taxdump_path3
    file x from refgenomes_ch1
    file x from refgenomesGB_ch1

    output:
    file 'Genomes.taxomonomy' into taxonomy_ch
    file '*.fna' into genomes_ch
    file 'Genome-downloader.log' into log_ch

    //script
    script:

    """
    #Delete false Genbak files
    rm -rf FALSE*
    #log part
    echo "Genome-downloader started at `date`" > Genome-downloader.log
    echo "Genomes dowloaded: " >> Genome-downloader.log
    find *.fna | wc -l >> Genome-downloader.log
    #copy part
    find *.fna | cut -f1 -d'-' > GC.list
    fetch-tax.pl GC.list  --taxdir=\$(<taxdump_path.txt) --item-type=taxid --levels=phylum class order family genus species
    mv GC.tax Genomes.taxomonomy
    for f in `cat GC.list`; do cp \$f-abbr.fna \$f-abbreviate.fna; done
    """
}
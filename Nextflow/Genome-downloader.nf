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
    --taxolevel      Taxonomic level of group - Choice between five taxa levels: phylum, class, order, family, genus, species          

    Optional arguments:
    --taxdump 		        Path to taxdump directory - automatic setup by default
    --refseq                Download refseq genomes, yes or no, activated by default
    --genbank               Download GenBank genomes, yes or no, activated by default
    --dRep                  dRep dereplication for group of interest, not activated by default - yes or no
    --ignoreGenomeQuality   don't run genome quality during dRep dereplication, yes or no, default = no
    --abbr                  Abreviate the deflines of the genomes, yes or no, default = no
    --prot                  Get proteins of corresponding genomes, yes or no, default = no
    --cpu                   number of cpus to use, default = 1
    
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
	exit 1, "Taxonomic levels: Choice between four taxa levels: phylum, class, order, family, genus, species "
}

//Path to companion
params.companion = '/opt/COMPANION/Genome-downloader_companion.py'

//Path to taxdump
params.taxdump = 'local'

//GenBank activated by default
params.refseq = 'yes'

//GenBank activated by default
params.genbank = 'yes'

//ignoreGenomeQuality
params.ignoreGenomeQuality = 'no'

//Dereplications
params.dRep = 'no'

//abbr
params.abbr = 'no'

//prot
params.prot = 'no'

//Not specified by user
//Path to project dir taxdump
taxdir = "$workflow.projectDir" + '/taxdump'
workingdir = file(taxdir)

//cpu
params.cpu = '1'

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
    file "ftp.sh" into refseq_ftp2
    file "GCF.tax" into refseq_tax1
    file 'Genome-downloader.log' into log_ch1

    //script
    script:
    if (params.refseq == 'yes'){
        println "Add Refseq Genomes activated"
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O refseq_sum.txt
        $companion refseq_sum.txt --mode=sum
        grep -v "#" refseq_sum-filt.txt | cut -f1 > GCF.list
        fetch-tax.pl GCF.list  --taxdir=\$(<taxdump_path.txt) --item-type=taxid --levels=phylum class order family genus species
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
        echo "RefSeq metadata" >> Genome-downloader.log
        """
    }
    else {
        println "Add Refseq Genomes NOT activated"
        """
        echo "Add Refseq Genomes NOT activated" > GCF.tax
        echo "Add Refseq Genomes NOT activated" > ftp.sh
        echo "Add Refseq metadata NOT activated" >> Genome-downloader.log
        """
    }
}

//Download Genbank metadata, compute taxonomy file and prudce download ftp file: OPTIONAL
process GenBank {

    //informations

	//input output
    input:
    file taxdump from taxdump_path2
    val companion from params.companion
    file 'Genome-downloader.log' from log_ch1
    
    output:
    file "GCA-ftp.sh" into genbank_ftp1
    file "GCA-ftp.sh" into genbank_ftp2
    file "GCA.tax" into genbank_tax1
    file 'Genome-downloader.log' into log_ch2

    //script
    script:
    if (params.genbank == 'yes'){
        println "Add GenBank Genomes activated"
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -O genbank_sum.txt
        $companion genbank_sum.txt --mode=sum
        grep -v "#" genbank_sum-filt.txt | cut -f1 > GCA.list
        fetch-tax.pl GCA.list  --taxdir=\$(<taxdump_path.txt) --item-type=taxid --levels=phylum class order family genus species
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
         echo "Add GenBank metadata activated" >> Genome-downloader.log
        """
    }
    else {
        println "Add GenBank Genomes NOT activated"
        """
        echo "Add GenBank Genomes NOT activated" > GCA.tax
        echo "Add GenBank Genomes NOT activated" > GCA-ftp.sh
        echo "Add GenBank metadata NOT activated" >> Genome-downloader.log
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
    file 'Genome-downloader.log' from log_ch2
    
    output:
    //file '*-abbr.fna' into refgenomes_ch1
    file '*.fna' into refgenomes_ch1
    file 'reduce-ftp.sh' into reduceRefFtp_ch
    file 'GCF.refgroup.uniq' into refgroupRefseqGC_ch
    file 'Genome-downloader.log' into log_ch3

    //script
    script:
    if (params.refseq == 'yes'){
        if (params.abbr == 'yes'){
            """
            #Produce list of GCF IDs with reference group and taxa levels
            $companion GCF.tax --mode=fetch --taxa=$taxa --refgroup=$group
            for f in `cat GCF.refgroup.uniq`; do grep \$f ftp.sh; done > reduce-ftp.sh
            bash reduce-ftp.sh
            find *.gz -type f -empty -print -delete
            gunzip *.gz
            find *.fna | cut -f1,2 -d"." > fna.list
            for f in `cat fna.list`; do inst-abbr-ids.pl \$f*.fna --id-regex=:DEF --id-prefix=\$f; done
            echo "Add RefSeq Genomes, abbr mode" >> Genome-downloader.log
            """
        }
        else if (params.abbr == 'no'){
            """
            #Produce list of GCF IDs with reference group and taxa levels
            $companion GCF.tax --mode=fetch --taxa=$taxa --refgroup=$group
            for f in `cat GCF.refgroup.uniq`; do grep \$f ftp.sh; done > reduce-ftp.sh
            bash reduce-ftp.sh
            find *.gz -type f -empty -print -delete
            gunzip *.gz
            echo "Add RefSeq Genomes, non abbr mode" >> Genome-downloader.log
            """
        }
    }
    else {
        println "Add RefSeq Genomes NOT activated"
        """
        echo "Add RefSeq Genomes NOT activated" > FALSER-abbr.fna
        echo "Add RefSeq Genomes NOT activated" > reduce-ftp.sh
        echo "GCF_FALSE" > GCF.refgroup.uniq
        echo "Add RefSeq Genomes NOT activated" >> Genome-downloader.log
        """
    }
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
    file 'Genome-downloader.log' from log_ch3
    
    output:
    file '*.fna' into refgenomesGB_ch1
    file 'GCA-reduce-ftp.sh' into reduceRefFtpGB_ch
    file 'Genome-downloader.log' into log_ch4

    //script
    script:
    if (params.genbank == 'yes'){
        if (params.abbr == 'yes'){
        println "Add GenBank Genomes activated"
            """
            #Produce list of GCA IDs with reference group and taxa levels
            $companion GCA.tax --mode=fetch --taxa=$taxa --refgroup=$group
            for f in `cat GCA.refgroup.uniq`; do grep \$f GCA-ftp.sh; done > GCA-reduce-ftp.sh
            bash GCA-reduce-ftp.sh
            find *.gz -type f -empty -print -delete
            gunzip *.gz
            find *.fna | cut -f1,2 -d"." > fna.list
            for f in `cat fna.list`; do inst-abbr-ids.pl \$f*.fna --id-regex=:DEF --id-prefix=\$f; done
            #for fix and proceed , false genbank files
            echo "Add GenBank Genomes activated" > FALSE-abbr.fna
            echo "Add GenBank Genomes activated" > FALSE-GCA-reduce-ftp.sh
            echo "Add GenBank Genomes activated, abbr mode" >> Genome-downloader.log
            """
        }
        else if (params.abbr == 'no'){
        println "Add GenBank Genomes activated"
            """
            #Produce list of GCA IDs with reference group and taxa levels
            $companion GCA.tax --mode=fetch --taxa=$taxa --refgroup=$group
            for f in `cat GCA.refgroup.uniq`; do grep \$f GCA-ftp.sh; done > GCA-reduce-ftp.sh
            bash GCA-reduce-ftp.sh
            find *.gz -type f -empty -print -delete
            gunzip *.gz
            #for fix and proceed , false genbank files
            echo "Add GenBank Genomes activated" > FALSE-abbr.fna
            echo "Add GenBank Genomes activated" > FALSE-GCA-reduce-ftp.sh
            echo "Add GenBank Genomes activated, non abbr mode" >> Genome-downloader.log
            """
        }
    }
    else {
        println "Add GenBank Genomes NOT activated"
        """
        echo "Add GenBank Genomes NOT activated" > FALSE-abbr.fna
        echo "Add GenBank Genomes NOT activated" > GCA-reduce-ftp.sh
        echo "Add GenBank Genomes NOT activated" >> Genome-downloader.log
        """
    }

}

//Dereplication with drep : OPTIONAL
process Dereplication {

    //informations

	//input output
    input:
    file x from refgenomes_ch1
    file x from refgenomesGB_ch1
    val cpu from params.cpu
    file 'Genome-downloader.log' from log_ch4
    
    output:
    file 'DEREPLICATED' into drep_ch1
    file 'DEREPLICATED' into drep_ch2
    file 'Genome-downloader.log' into log_ch
    file 'Genome-downloader.log' into log_ch5

    //script
    script:
    if (params.dRep == 'yes'){
        println "Genomes dereplication activated"
        if (params.ignoreGenomeQuality == 'no') {
            """
            #Delete false Genbak files
            rm -rf FALSE*
            mkdir GEN
            mv *.fna GEN/
            dRep dereplicate DREP -g GEN/*.fna -p $cpu
            mkdir DEREPLICATED/
            mv DREP/dereplicated_genomes/*.fna DEREPLICATED/
            echo "Drep Dereplication activated" >> Genome-downloader.log
            """
        }
        else if (params.ignoreGenomeQuality == 'yes') {
            """
            #Delete false Genbak files
            rm -rf FALSE*
            mkdir GEN
            mv *.fna GEN/
            dRep dereplicate DREP -g GEN/*.fna -p $cpu --ignoreGenomeQuality
            mkdir DEREPLICATED/
            mv DREP/dereplicated_genomes/*.fna DEREPLICATED/
            echo "Drep Dereplication activated" >> Genome-downloader.log
            """
        }
    }
    else {
        println "Genomes dereplication not activated"
        """
        #Delete false Genbak files
        rm -rf FALSE*
        mkdir DEREPLICATED/
        mv *.fna DEREPLICATED/
        echo "Drep Dereplication not activated" >> Genome-downloader.log
        """

    }

}

//Get Prot : OPTIONAL
process GetProt {

    //informations

	//input output
    input:
    file 'DEREPLICATED' from drep_ch2
    file "ftp.sh" from refseq_ftp2
    file "GCA-ftp.sh" from genbank_ftp2
    file 'Genome-downloader.log' from log_ch5
    
    output:
    file 'PROT' into prot_ch
    file 'Genome-downloader.log' into log_ch6

    //script
    script:
    if (params.prot == 'yes'){
        println "Prot download activated"
        """
        #log part
        echo test > Genome-downloader.log
        #Merge ftp
        cat ftp.sh GCA-ftp.sh > merge-ftp.sh
        grep -v 'activated' merge-ftp.sh > temp; mv temp merge-ftp.sh
        #Collecte GCA/F number
        cp DEREPLICATED/*.fna .
        find *.fna > fna.list
        sed -i -e 's/-abbr.fna//g' fna.list
        sed -i -e 's/.fna//g' fna.list
        rm -rf *.fna
        #Get ftp part
        for f in `cat fna.list`; do grep \$f merge-ftp.sh; done > prot-ftp.sh
        sed -i -e 's/_genomic.fna.gz/_protein.faa.gz/g' prot-ftp.sh
        sed -i -e 's/fna.gz/faa.gz/g' prot-ftp.sh
        bash prot-ftp.sh
        find *.gz -type f -empty -print -delete
        gunzip *.gz
        mkdir PROT
        mv *.faa PROT/
        """
    }
    else {
        println "Prot download NOT activated"
        """
        #log part
        echo "Prot download NOT activated" > Genome-downloader.log
        mkdir PROT
        echo "Prot download NOT activated" > PROT/info.faa
        """

    }

}


process CombineGenomes {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file taxdump from taxdump_path3
    file 'DEREPLICATED' from drep_ch1
    file 'PROT' from prot_ch
    file 'Genome-downloader.log' from log_ch6

    output:
    file 'Genomes.taxomonomy' into taxonomy_ch
    file 'GENOMES' into genomes_ch
    file 'PROTEINS' into proteins_ch
    file 'Genome-downloader.log' into logfinal_ch

    //script
    script:

    """
    mv DEREPLICATED/*.fna .
    #log part
    echo "Genome-downloader started at `date`" > Genome-downloader.log
    echo "Genomes dowloaded: " >> Genome-downloader.log
    find *.fna | wc -l >> Genome-downloader.log
    echo "Genome-downloader, version: " >> Genome-downloader.log
    echo "3.0.0 " >> Genome-downloader.log
    #copy part
    #find *.fna | cut -f1 -d'-' > GC.list
    find *.fna  > GC.list
    mkdir GENOMES/
    mv *.fna GENOMES/
    sed -i -e 's/.fna//g' GC.list
    fetch-tax.pl GC.list  --taxdir=\$(<taxdump_path.txt) --item-type=taxid --levels=phylum class order family genus species
    mv GC.tax Genomes.taxomonomy
    mkdir PROTEINS
    mv PROT/*.faa PROTEINS/
    """
}

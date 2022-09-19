#!/usr/bin/env nextflow
/*
========================================================================================
                         Nextflow-ORganismPlacER(ORPER)
========================================================================================
GIT url : https://github.com/Lcornet/ORPER
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
    ORPER performs a phylogenetic placement of SSU sequences in a tree, composed of RefSeq genomes, and 
    constrained by a ribosomal phylogenomic tree.
    
    Citation:
    Please cite : TODO

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow run ORPER-durandal.nf --GenBank=yes --refgroup=Nostocales --outgroup=Gloeobacterales --taxa=order 
    --ribodb=PATH-to-ribodb-files/ --companion=/opt/ORPER-companion.py --cpu=20 --SSU=Nostocales-silva.fasta 
    --cdhit=yes --drep=yes
    
    Mandatory arguments:
    --refgroup          RefSeq group of interest
    --outgroup          Refseq outgroup
    --reftaxolevel      Taxonomic level of reference group - Choice between four taxa levels: phylum, class, order, family
    --outtaxolevel      Taxonomic level of out group - Choice between four taxa levels: phylum, class, order, family
    --SSU               Path to fasta file containing SSU sequences                     

    Optional arguments:
    --ribodb            Path to directory containing ribodb fasta files
    --companion         Path to ORPER-companion.py 
    --taxdump 		    Path to taxdump directory, created if not set
    --cpu               Number of cpus, default = 1
    --refgenbank        add GenBank for reference group, Deactivated by default - yes or no
    --outgenbank        add GenBank for out group, Deactivated by default - yes or no
    --outgenbank        Download GenBank metadata, activated by default - yes or no. Needed if refgenbank or outgenbank is activated
    --dRep              dRep dereplication of RefSeq group of interest, activated by default - yes or no
    --cdhit             cdhit dereplication of provided SSU sequences, activated by default with 100% threshold - yes or no
    --shrink            TreeShrink cutoff value
    
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
params.refgroup = null
if (params.refgroup == null) {
	exit 1, "RefSeq group mandatory,RefSeq group used to compute backdone of the tree"
}

//Name of the Outgroup : Mandatory
params.outgroup = null
if (params.outgroup == null) {
	exit 1, "RefSeq group mandatory,RefSeq group used for outgroup"
}

//Path to SSU fasta file : Mandatory
params.SSU = null
if (params.SSU == null) {
	exit 1, "Path to fasta file containing SSU sequences"
}

//Choice between four taxa levels for refgroup: phylum, class, order, family: Mandatory
params.reftaxolevel = null
if (params.reftaxolevel == null) {
	exit 1, "Path to fasta file containing SSU sequences"
}

//Choice between four taxa levels for outgroup: phylum, class, order, family: Mandatory
params.outtaxolevel = null
if (params.outtaxolevel == null) {
	exit 1, "Path to fasta file containing SSU sequences"
}

//Path to companion
params.companion = '/opt/ORPER/ORPER-companion.py'

//Ribodb
params.ribodb = 'local'

//Path to taxdump
params.taxdump = 'local'

//Number of cpus
params.cpu = '1'

//GenBank metadata activated by default
params.genbank = 'yes'

//GenBank genomes for refgroup deactivated by default
params.refgenbank = 'no'

//GenBank genomes for outgroup deactivated by default
params.outgenbank = 'no'

//Treeshrink value
params.shrink = '0.1'

//Not specified by user
//Path to project dir taxdump
taxdir = "$workflow.projectDir" + '/taxdump'
workingdir = file(taxdir)

//Path to project dir RiboDB
ribodir = "$workflow.projectDir" + '/ORPER-ribodb'
workingribo = file(ribodir)

//outdir
params.outdir='ORPER-results'

//Dereplications
params.dRep = 'yes'
params.cdhit = 'yes'


/*
CORE PROGRAM
*/

//Load input files
taxdump_ch = Channel.fromPath(params.taxdump)
ribodb_ch = Channel.fromPath(params.ribodb)
ssu_ch = Channel.fromPath(params.SSU)

//Get ribodb files
process RiboddbSetUp {
	//informations

	//input output
    input:
    val ribodb from ribodb_ch
    
    output:
    val ribodir into ribodir_ch1
    

    //script
    script:

    ribo = 'na'

    if (params.ribodb == 'local'){
        println "ORPER-INFO: RiboDB not specified -> project dir"

        if( !workingribo.exists() ) {
            println "ORPER-INFO: RiboDB dir not found in project dir -> Created"
            if( !workingribo.mkdirs() )    {
                exit 1, "Cannot create working directory"
            }

            ribodir = workingribo

            """
            cd $workingribo
            git clone https://bitbucket.org/phylogeno/42-ribo-msas/src/master/MSAs/prokaryotes/
            mv prokaryotes/MSAs/prokaryotes/*.ali .
            ali2fasta.pl *.ali --degap 2> log
            find *.fasta | cut -f1 -d'.' > ribo.list
            #for f in `cat ribo.list`; do mv \$f.fasta \$f-prot_abbr.fasta; done
            rm -f *.ali
            echo $workingribo > ribodb_path.txt
            """
        }
        else {
            println "ORPER-INFO: Taxdump dir found in project dir -> Used"

            ribodir = workingribo 

 	        """
            echo $workingribo > ribodb_path.txt
		    """           

        }
    }

	else{
        println "ORPER-INFO: RiboDB specified"

        ribodir = ribodb

		"""
        echo $ribodb > ribodb_path.txt
		"""		
    }
}


//Taxonomy, set taxdump if not specifed
process Taxonomy {
	//informations

	//input output
    input:
    val taxdump from taxdump_ch
    
    output:
    file "taxdump_path.txt" into taxdump_path1
    file "taxdump_path.txt" into taxdump_path2
    val taxdir into taxdir_ch1
    val taxdir into taxdir_ch2

    //script
    script:

    taxdir = 'na'

    if (params.taxdump == 'local'){
        println "ORPER-INFO: Taxdump not specified -> project dir"

        if( !workingdir.exists() ) {
            println "ORPER-INFO: Taxdump dir not found in project dir -> Created"
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
            println "ORPER-INFO: Taxdump dir found in project dir -> Used"

            taxdir = workingdir 

 	        """
            echo $workingdir > taxdump_path.txt
		    """           

        }
    }

	else{
        println "ORPER-INFO: Taxdump specified"

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
    file "GCF.tax" into refseq_tax2

    //script
    script:

    """
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O refseq_sum.txt
    $companion refseq_sum.txt --mode=sum
    grep -v "#" refseq_sum-filt.txt | cut -f1 > GCF.list
    fetch-tax.pl GCF.list  --taxdir=\$(<taxdump_path.txt) --item-type=taxid --levels=phylum class order family
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
    file "GCA-ftp.sh" into genbank_ftp2
    file "GCA.tax" into genbank_tax1
    file "GCA.tax" into genbank_tax2

    //script
    script:
    if (params.genbank == 'yes'){
        println "Add GenBank Genomes activated"
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -O genbank_sum.txt
        $companion genbank_sum.txt --mode=sum
        grep -v "#" genbank_sum-filt.txt | cut -f1 > GCA.list
        fetch-tax.pl GCA.list  --taxdir=\$(<taxdump_path.txt) --item-type=taxid --levels=phylum class order family
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
process GetRefGenomesRefseq {
	//informations

	//input output
    input:
    val refgroup from params.refgroup
    val taxa from params.reftaxolevel
    val companion from params.companion
    file "ftp.sh" from refseq_ftp1
    file "GCF.tax" from refseq_tax1
    
    output:
    file '*-abbr.fna' into refgenomes_ch1
    file '*-abbr.fna' into refgenomes_ch2
    file '*-abbr.fna' into refgenomes_ch3
    file 'reduce-ftp.sh' into reduceRefFtp_ch
    file 'GCF.refgroup.uniq' into refgroupRefseqGC_ch

    //script
    script:

    """
    #Produce list of GCF IDs with reference group and taxa levels
    $companion GCF.tax --mode=fetch --taxa=$taxa --refgroup=$refgroup
    for f in `cat GCF.refgroup.uniq`; do grep \$f ftp.sh; done > reduce-ftp.sh
    bash reduce-ftp.sh
    gunzip *.gz
    find *.fna | cut -f1,2 -d"." > fna.list
    for f in `cat fna.list`; do inst-abbr-ids.pl \$f*.fna --id-regex=:DEF --id-prefix=\$f; done
    """
}

//Get GenBank genome for the reference group, abbr files: OPTIONAL
process GetRefGenomesGenbank {

    //informations

	//input output
    input:
    val refgroup from params.refgroup
    val taxa from params.reftaxolevel
    val companion from params.companion
    file "GCA-ftp.sh" from genbank_ftp1
    file "GCA.tax" from genbank_tax1
    file 'GCF.refgroup.uniq' from refgroupRefseqGC_ch
    
    output:
    file '*-abbr.fna' into refgenomesGB_ch1
    file '*-abbr.fna' into refgenomesGB_ch2
    file '*-abbr.fna' into refgenomesGB_ch3
    file 'GCA-reduce-ftp.sh' into reduceRefFtpGB_ch

    //script
    script:
    if (params.refgenbank == 'yes'){
        println "Add GenBank Genomes activated"
        """
        #Produce list of GCA IDs with reference group and taxa levels
        $companion GCA.tax --mode=fetch --taxa=$taxa --refgroup=$refgroup
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

//run contamination evaluation
process RefGenomesCheckm {
	//informations

	//input output
    input:
    file x from refgenomes_ch1
    file x from refgenomesGB_ch1
    val cpu from params.cpu

    output:
    file "RefGenomes.Checkm" into refGenomesCheckm_ch

    //script
    script:

    """
    #Delete false Genbak files
    rm -rf FALSE*
    mkdir RefGenomes
    mv *.fna RefGenomes/
    checkm lineage_wf -t $cpu -x fna RefGenomes runc > checkm.result
    echo "#genome,completeness,contamination" > part1
    tr -s " " < checkm.result | grep "GC" | cut -f2,14,15 -d" " > part2
    sed -i -e 's/ /,/g' part2
    cat part1 part2 > RefGenomes.Checkm
    """
}

//Run rnammer on all genomes
process RefGenomesBarnapp  {
	//informations

	//input output
    input:
    file x from refgenomes_ch2
    file x from refgenomesGB_ch2
    val companion from params.companion

    output:
    file "genome-with-ssu.list" into refGenomeWithSsu_ch
    file "all_16s-nodupe.fna" into refRnammer_ch

    //script
    script:

    """
    #Delete false Genbak files
    rm -rf FALSE*
    find *.fna | cut -f1 -d"-" > fna.list
    for f in `cat fna.list`; do barrnap \$f-abbr.fna --outseq \$f-barnap.fna --threads 1; done
    #for f in `cat fna.list`; do fasta2ali.pl \$f-barnap.fna; done
    #for f in `cat fna.list`; do grep -A1 '16S' \$f-barnap.ali > \$f-16s.ali; done
    #for f in `cat fna.list`; do ali2fasta.pl \$f-16s.ali; mv \$f-16s.fasta \$f-16s.fna; done
    cat *barnap.fna > all_barnap.fna
    $companion all_barnap.fna --mode=barnap
    #cat *-16s.fna > all_16s.fna
    #$companion all_16s.fna --mode=barnap
    #RNAMMER
    #for f in `cat fna.list`; do rnammer -S bac -m ssu -d -gff \$f.gff -h \$f.hmm -f `basename \$f .fa`_16s.fna < \$f-abbr.fna; done
    #cat *_16s.fna > all_16s.fna
    #$companion all_16s.fna --mode=rnammer
    """
}


//Filter genomes based on completeness, contamination and 16s presence
process RefGenomesFilter {
	//informations

	//input output
    input:
    file "genome-with-ssu.list" from refGenomeWithSsu_ch
    file "RefGenomes.Checkm" from refGenomesCheckm_ch
    val companion from params.companion

    output:
    file "reliable-genomes.list" into refReliablegenomes_ch

    //script
    script:

    """
    $companion RefGenomes.Checkm --mode=checkm --ssu=yes
    """
}

//Dereplication with drep : OPTIONAL
process RefGenomesDereplication {

    //informations

	//input output
    input:
    file x from refgenomes_ch3
    file x from refgenomesGB_ch3
    file "reliable-genomes.list" from refReliablegenomes_ch
    val cpu from params.cpu
    
    output:
    file "reliable-genomes-dereplicated.list" into refDrepReliablegenomes_ch

    //script
    script:
    if (params.dRep == 'yes'){
        println "Ref Genomes dereplication activated"
        """
        mkdir Genomes
        for f in `cat reliable-genomes.list`; do mv \$f*.fna Genomes; done
        rm -rf *.fna
        dRep dereplicate DREP -g Genomes/*.fna -p $cpu
        cd DREP/dereplicated_genomes/
        find *.fna | cut -f1 -d'-' > reliable-genomes-dereplicated.list
        mv reliable-genomes-dereplicated.list ../../
        """
    }
    else {
        println "Ref Genomes dereplication not activated"
        """
        cp reliable-genomes.list reliable-genomes-dereplicated.list
        """

    }

}

//Get proteomes of reliable genomes
process GetRefRelProteomes {
	//informations

	//input output
    input:
    file "reliable-genomes.list" from refDrepReliablegenomes_ch
    file 'reduce-ftp.sh' from reduceRefFtp_ch
    file 'GCA-reduce-ftp.sh' from reduceRefFtpGB_ch

    output:
    file '*-abbr.faa' into refReliableproteomes_ch

    //script
    script:

    """
    cat reduce-ftp.sh GCA-reduce-ftp.sh > combined-reduce-ftp.sh
    for f in `cat reliable-genomes.list `; do grep "\$f" combined-reduce-ftp.sh ; done > reliable.sh
    sed -i -e 's/_genomic.fna/_protein.faa/g' reliable.sh
    sed -i -e 's/.fna.gz/.faa.gz/g' reliable.sh
    bash reliable.sh
    find . -name '*.gz' -size 0 | cut -f2 -d"/" > empty.list
    for f in `cat empty.list `; do rm -rf \$f; done
    gunzip *.gz
    find *.faa | cut -f1,2 -d"." > faa.list
    for f in `cat faa.list`; do inst-abbr-ids.pl \$f*.faa --id-regex=:DEF --id-prefix=\$f; done
    """
}

//OutGroup part
//Get Refseq genome for the out group, abbr files
process GetOutGenomesRefSeq {
	//informations

	//input output
    input:
    val outgroup from params.outgroup
    val taxa from params.outtaxolevel
    val companion from params.companion
    file "ftp.sh" from refseq_ftp2
    file "GCF.tax" from refseq_tax2
    
    output:
    file '*-abbr.fna' into outgenomes_ch1
    file '*-abbr.fna' into outgenomes_ch2
    file '*-abbr.fna' into outgenomes_ch3
    file 'reduce-ftp.sh' into reduceOutFtp_ch
    file 'GCF.outgroup.uniq' into outgroupRefseqGC_ch

    //script
    script:

    """
    #Produce list of GCF IDs with reference outgroup and taxa levels
    $companion GCF.tax --mode=fetch --taxa=$taxa --refgroup=$outgroup
    mv GCF.refgroup.uniq GCF.outgroup.uniq
    for f in `cat GCF.outgroup.uniq`; do grep \$f ftp.sh; done > reduce-ftp.sh
    bash reduce-ftp.sh
    #create a false file, gunziped, in case no outgroup could be found in refseq
    echo ">FALSE FALSE" > False.1.fna
    gzip False.1.fna
    #all files together, including false
    gunzip *.gz
    find *.fna | cut -f1,2 -d"." > fna.list
    for f in `cat fna.list`; do inst-abbr-ids.pl \$f*.fna --id-regex=:DEF --id-prefix=\$f; done
    """
}

//Get Genbank genome for the out group, abbr files: Optional
process GetOutGenomesGenbank {

    //informations

	//input output
    input:
    val outgroup from params.outgroup
    val taxa from params.outtaxolevel
    val companion from params.companion
    file "GCA-ftp.sh" from genbank_ftp2
    file "GCA.tax" from genbank_tax2
    file 'GCF.outgroup.uniq' from outgroupRefseqGC_ch
    
    output:
    file '*-abbr.fna' into outgenomesGB_ch1
    file '*-abbr.fna' into outgenomesGB_ch2
    file '*-abbr.fna' into outgenomesGB_ch3
    file 'GCA-reduce-ftp.sh' into reduceOutFtpGB_ch

    //script
    script:
    if (params.outgenbank == 'yes'){
        println "Add GenBank Genomes activated"
        """
        #Produce list of GCA IDs with outgroup group and taxa levels
        cp GCF.outgroup.uniq GCF.refgroup.uniq #for companion
        $companion GCA.tax --mode=fetch --taxa=$taxa --refgroup=$outgroup
        mv GCA.refgroup.uniq GCA.outgroup.uniq
        for f in `cat GCA.outgroup.uniq`; do grep \$f GCA-ftp.sh; done > GCA-reduce-ftp.sh
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

//run contamination evaluation
process OutGenomesCheckm {
	//informations

	//input output
    input:
    file x from outgenomes_ch1
    file x from outgenomesGB_ch1
    val cpu from params.cpu

    output:
    file "OutGenomes.Checkm" into outGenomesCheckm_ch

    //script
    script:

    """
    #pyenv local 2.7.6
    #delete false file
    rm -f False*
    rm -f FALSE*
    mkdir OutGenomes
    mv *.fna OutGenomes/
    checkm lineage_wf -t $cpu -x fna OutGenomes runc > checkm.result
    echo "#genome,completeness,contamination" > part1
    tr -s " " < checkm.result | grep "GC" | cut -f2,14,15 -d" " > part2
    sed -i -e 's/ /,/g' part2
    cat part1 part2 > OutGenomes.Checkm
    """
}

//Run rnammer on all genomes
process OutGenomesBarnap  {
	//informations

	//input output
    input:
    file x from outgenomes_ch2
    file x from outgenomesGB_ch2
    val companion from params.companion

    output:
    file "genome-with-ssu.list" into outGenomeWithSsu_ch
    file "all_16s-nodupe.fna" into outRnammer_ch

    //script
    script:

    """
    #delete false file
    rm -f False*
    rm -f FALSE*
    find *.fna | cut -f1 -d"-" > fna.list
    for f in `cat fna.list`; do barrnap \$f-abbr.fna --outseq \$f-barnap.fna --threads 1; done
    #for f in `cat fna.list`; do fasta2ali.pl \$f-barnap.fna; done
    #for f in `cat fna.list`; do grep -A1 '16S' \$f-barnap.ali > \$f-16s.ali; done
    #for f in `cat fna.list`; do ali2fasta.pl \$f-16s.ali; mv \$f-16s.fasta \$f-16s.fna; done
    #cat *-16s.fna > all_16s.fna
    cat *barnap.fna > all_barnap.fna
    $companion all_barnap.fna --mode=barnap
    """
}

//Filter genomes based on completeness, contamination and 16s presence
process OutGenomesFilter {
	//informations

	//input output
    input:
    file "genome-with-ssu.list" from outGenomeWithSsu_ch
    file "OutGenomes.Checkm" from outGenomesCheckm_ch
    val companion from params.companion

    output:
    file "reliable-genomes.list" into outReliablegenomes_ch

    //script
    script:

    """
    $companion OutGenomes.Checkm --mode=checkm --ssu=yes
    """
}

//Get proteomes of reliable genomes
process GetOutRelProteomes {
	//informations

	//input output
    input:
    file "reliable-genomes.list" from outReliablegenomes_ch
    file 'reduce-ftp.sh' from reduceOutFtp_ch
    file 'GCA-reduce-ftp.sh' from reduceOutFtpGB_ch
    file x from outgenomes_ch3
    file x from outgenomesGB_ch3

    output:
    file '*-abbr.faa' into outReliableproteomes_ch

    //script
    script:

    """
    #Take only 10 Genomes for outgroup
    head -n10  reliable-genomes.list > reliable-n10.list
    mkdir N10
    for f in `cat reliable-n10.list`; do mv \$f*.fna N10/; done
    rm -rf *.fna
    cd N10/
    for f in `cat ../reliable-n10.list`; do prodigal -i \$f-abbr.fna -o \$f.genes -a \$f.faa; done
    for f in `cat ../reliable-n10.list`; do inst-abbr-ids.pl \$f.faa --id-regex=:DEF; done
    for f in `cat ../reliable-n10.list`; do sed -i -e 's/|GCA_/GCA_/g' \$f-abbr.faa; done
    for f in `cat ../reliable-n10.list`; do sed -i -e 's/|GCF_/GCF_/g' \$f-abbr.faa; done
    mv *abbr.faa ../
    cd  ../
    """
}

//Common part :  Ref and Out group
//Run forty to enrich ribodb with reference genomes
process RiboDBFortytwo {
	//informations

	//input output
    input:
    val cpu from params.cpu
    file x from refReliableproteomes_ch
    file x from outReliableproteomes_ch
    val ribodb from params.ribodb
    val taxdir from taxdir_ch1
    val companion from params.companion
    val ribodir from ribodir_ch1

    output:
    file '*enrich.fasta' into enrichedRiboDB_ch

    //script
    script:

    """
    #Reference organism part 
    mkdir ref-banks
    cd ref-banks/
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/022/565/GCA_000022565.1_ASM2256v1/GCA_000022565.1_ASM2256v1_protein.faa.gz -O GCA_000022565.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/019/605/GCA_000019605.1_ASM1960v1/GCA_000019605.1_ASM1960v1_protein.faa.gz -O GCA_000019605.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/956/175/GCA_000956175.1_ASM95617v1/GCA_000956175.1_ASM95617v1_protein.faa.gz -O GCA_000956175.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/024/005/GCA_000024005.1_ASM2400v1/GCA_000024005.1_ASM2400v1_protein.faa.gz -O GCA_000024005.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/905/GCA_000026905.1_ASM2690v1/GCA_000026905.1_ASM2690v1_protein.faa.gz -O GCA_000026905.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/545/GCA_000026545.1_ASM2654v1/GCA_000026545.1_ASM2654v1_protein.faa.gz -O GCA_000026545.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/145/985/GCA_000145985.1_ASM14598v1/GCA_000145985.1_ASM14598v1_protein.faa.gz -O GCA_000145985.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/165/505/GCA_000165505.1_ASM16550v1/GCA_000165505.1_ASM16550v1_protein.faa.gz -O GCA_000165505.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/085/GCA_000008085.1_ASM808v1/GCA_000008085.1_ASM808v1_protein.faa.gz -O GCA_000008085.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/698/785/GCA_000698785.1_ASM69878v1/GCA_000698785.1_ASM69878v1_protein.faa.gz -O GCA_000698785.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/725/425/GCA_000725425.1_ASM72542v1/GCA_000725425.1_ASM72542v1_protein.faa.gz -O GCA_000725425.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/725/GCA_000091725.1_ASM9172v1/GCA_000091725.1_ASM9172v1_protein.faa.gz -O GCA_000091725.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/011/505/GCA_000011505.1_ASM1150v1/GCA_000011505.1_ASM1150v1_protein.faa.gz -O GCA_000011505.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/012/285/GCA_000012285.1_ASM1228v1/GCA_000012285.1_ASM1228v1_protein.faa.gz -O GCA_000012285.1.faa.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/585/GCA_000014585.1_ASM1458v1/GCA_000014585.1_ASM1458v1_protein.faa.gz -O GCA_000014585.1.faa.gz
    gunzip *.gz
    find *.faa | cut -f1,2 -d"." > faa.list
    for f in `cat faa.list`; do makeblastdb -in \$f.faa -dbtype prot -parse_seqids -out \$f; done
    for f in `cat faa.list `; do echo ".psq" ; done > end.list
    paste faa.list end.list > part1
    sed -i -e 's/\t//g' part1
    #paste part1 faa.list > ref-bank-mapper.idm
    paste faa.list part1 > ref-bank-mapper.idm
    cd ..

    #Define Queries
    echo "Sulfolobus solfataricus_2287" >> queries.idl
    echo "Thermoproteus uzoniensis_999630" >> queries.idl
    echo "Vulcanisaeta moutnovskia_985053" >> queries.idl
    echo "Flavobacterium psychrophilum_96345" >> queries.idl
    echo "Brucella suis_645170" >> queries.idl
    echo "Burkholderia mallei_13373" >> queries.idl
    echo "Neisseria meningitidis_487" >> queries.idl
    echo "Helicobacter pylori_210" >> queries.idl
    echo "Escherichia coli_83333" >> queries.idl
    echo "Yersinia pestis_632" >> queries.idl
    echo "Pseudomonas aeruginosa_287" >> queries.idl
    echo "Francisella philomiragia_28110" >> queries.idl
    echo "Xanthomonas citri_611301" >> queries.idl
    echo "Chlamydia pneumoniae_83558" >> queries.idl
    echo "Corynebacterium pseudotuberculosis_1719" >> queries.idl
    echo "Mycobacterium tuberculosis_1773" >> queries.idl
    echo "Bacillus anthracis_1392" >> queries.idl
    echo "Listeria monocytogenes_1639" >> queries.idl
    echo "Staphylococcus aureus_1074919" >> queries.idl
    echo "Streptococcus agalactiae_1311" >> queries.idl

    #Part for genomes to add
    mkdir genomes-to-add
    mv *abbr.faa genomes-to-add/
    cd genomes-to-add/
    find *.faa | cut -f1 -d"-" > genomes.list
    for f in `cat genomes.list`; do makeblastdb -in \$f-abbr.faa -dbtype prot -parse_seqids -out \$f; done
    for f in `cat genomes.list`; do echo ".psq"; done  > end.list
    paste genomes.list end.list > part1
    sed -i -e 's/\t//g' part1
    find *abbr.faa | cut -f1,2 -d"." > part2
    #paste part1 part2 > bank-mapper.idm
    paste part2 part1 > bank-mapper.idm
    cd ..
   
    #Part for alignements
    mkdir ribodb
    cp $ribodir/*.fasta ribodb/
    
    #Generate yaml
    yaml-generator-42.pl --run_mode=phylogenomic --out_suffix=-ORPER --queries queries.idl --evalue=1e-05 --homologues_seg=yes \
    --max_target_seqs=10000 --templates_seg=no --bank_dir genomes-to-add --bank_suffix=.psq --bank_mapper genomes-to-add/bank-mapper.idm --ref_brh=on \
    --ref_bank_dir ref-banks --ref_bank_suffix=.psq --ref_bank_mapper ref-banks/ref-bank-mapper.idm --ref_org_mul=0.3 --ref_score_mul=0.99 \
    --trim_homologues==off --ali_keep_lengthened_seqs=keep --aligner_mode=off --tax_reports=off --tax_dir $taxdir \
    --megan_like --tol_check=off

    #run forty-two
    forty-two.pl ribodb/*.fasta --config=config-ORPER.yaml --verbosity=1 --threads=$cpu

    #extract new sequences from enriched ribodb
    cd ribodb/
    fasta2ali.pl *-ORPER.fasta
    grep -c "#NEW" *ORPER.fasta > count-enrich.list
    $companion count-enrich.list --mode=forty
    for f in `cat enrich.list`; do grep -A1 "#NEW#" \$f.ali > \$f-enrich.ali; done
    ali2fasta.pl --degap --noguessing *enrich.ali
    mv *enrich.fasta ../
    cd ..
    """
}

process AlignmentMafft {
	//informations

	//input output
    input:
    file '*enrich.fasta' from enrichedRiboDB_ch

    output:
    file '*enrich-ali.fasta' into alignments_ch

    //script
    script:

    """
    #delete empty file
    #find *enrich.fasta -size  0 -print -delete
    find *.fasta | cut -f1 -d"." > enrich.list
    #for f in `cat enrich.list`; do mafft --anysymbol --auto --reorder \$f.fasta > \$f-ali.fasta; done
    for f in `cat enrich.list`; do muscle3.8.31_i86linux64 -in \$f.fasta -out \$f-ali.fasta; done
    """
}


process ConcatScafos {
	//informations

	//input output
    input:
    file '*enrich-ali.fasta' from alignments_ch

    output:
    file 'data-ass.ali' into matrix_ch
    
    //script
    script:

    """
    #transform aligned fasta into ali format
    mkdir aligned
    mkdir a2p
    mv *.fasta aligned/
    cd aligned/
    fasta2ali.pl *.fasta
    #Filter with BMGE
    ali2phylip.pl *.ali --bmge-mask=medium --min=0.5 --ali
    mv *a2p* ../a2p/
    cd ../

    #Produce concat with scafos
    scafos in=a2p out=otu
    scafos in=a2p out=data otu=otu/otu-freq.otu o=ov
    scafos in=data out=data-ass otu=otu/otu-freq.otu gamma=yes o=gclv g=30 format=fpm
    cp data-ass/data-ass.ali .

    """
}

process ReferenceTreeRaxml {
	//informations

	//input output
    input:
    file 'data-ass.ali' from matrix_ch
    val cpu from params.cpu

    output:
    file 'reference.tre' into referenceTree_ch
    file 'GC.list' into referenceTreeList_ch

    //script
    script:

    """
    #transform matrix into phylip file  map-ids
    ali2phylip.pl data-ass.ali --map-ids
    #compute tree
    # -x -f a -N 100 will do a 100x rapid bootstrap analysis
    raxmlHPC-PTHREADS-AVX -T $cpu -s data-ass.phy -n data-ass-RAXML-PROTGAMMALGF-100xRAPIDBP -m PROTGAMMALGF -N 100 -f a -x 1975021703574 -p 1975021703574
    cp data-ass.idm RAxML_bipartitions.idm
    sed -i -e 's/-abbr//g' RAxML_bipartitions.idm
    cut -f1 -d"@" RAxML_bipartitions.idm > f1
    cut -f2 RAxML_bipartitions.idm > f2
    paste f1 f2 > RAxML_bipartitions.idm
    format-tree.pl RAxML_bipartitions.data-ass-RAXML-PROTGAMMALGF-100xRAPIDBP --map-ids
    mv RAxML_bipartitions.tre reference.tre
    tree2list.pl reference.tre
    grep -v "#" reference.idl > GC.list
    """
}

//Dereplication of provided 16s file : OPTIONAL
process SSUDereplication {

    //informations

	//input output
    input:
    file 'SSU.fasta' from ssu_ch
    val cpu from params.cpu
    
    output:
    file 'derepliacted-SSU.fasta' into ssuDereplicated_ch

    //script
    script:
    if (params.cdhit == 'yes'){
        println "SSU CD-HIT dereplication activated"
        """
        cd-hit -i SSU.fasta -o derepliacted-SSU.fasta -c 1.0 -T $cpu
        """
    }
    else {
        println "SSU CD-HIT dereplication not activated"
        """
        cp SSU.fasta derepliacted-SSU.fasta
        """

    }

}


//compute the 16s alignement and the constrained tree
process ConstrainTreeRaxml {
	//informations

	//input output
    input:
    file 'reference.tre' from referenceTree_ch
    file "Refall_16s-nodupe.fna" from refRnammer_ch
    file "Outall_16s-nodupe.fna" from outRnammer_ch
    file 'SSU.fasta' from ssuDereplicated_ch
    file 'GC.list' from referenceTreeList_ch
    val taxdir from taxdir_ch2
    val cpu from params.cpu
    val shrink from params.shrink
    val companion from params.companion

    output:
    file 'Constained-tree.nex' into constTreeNexus_ch
    file 'Constained-tree.tre' into constTreeTre_ch
    file 'SSU-combined-ali-a2p.fasta' into alignementBmge_ch
    file 'SSU-combined-ali.fasta' into alignement_ch
    file 'GC.list' into gcfList_ch


    //script
    script:

    """
    #Process input SSU fasta file, corrected unix fasta file
    fasta2ali.pl SSU.fasta
    ali2fasta.pl SSU.ali

    #Shorten sequence to first blank
    inst-abbr-ids.pl SSU.fasta --id-regex=:DEF
    sed -i -e 's/|//g' SSU-abbr.fasta

    #Process refgroup 16s sequences -- to replace in companion --
    #ensure that rejected sequence in scafos are not present in reference 16s file
    cat Refall_16s-nodupe.fna Outall_16s-nodupe.fna > all_16s-nodupe.fasta
    $companion all_16s-nodupe.fasta --mode=ConstrainTreeRaxml
    
    #combine reference ans SSU, fasta files
    cat all_16s-nodupe-list.fasta SSU-abbr.fasta > SSU-combined.fasta

    #align comibined file, use BMGE
    #mafft --adjustdirection --anysymbol --auto --reorder  SSU-combined.fasta > SSU-combined-ali.fasta
    muscle3.8.31_i86linux64 -in SSU-combined.fasta -out SSU-combined-ali.fasta
    fasta2ali.pl SSU-combined-ali.fasta
    sed -i -e 's/ //g' SSU-combined-ali.ali
    ali2phylip.pl SSU-combined-ali.ali --bmge-mask=medium --max=0.6 --ali
    ali2phylip.pl SSU-combined-ali-a2p.ali --p80
    ali2fasta.pl SSU-combined-ali-a2p.ali
    raxmlHPC-PTHREADS-AVX -T $cpu -r reference.tre -s SSU-combined-ali-a2p.p80 \
    -n SSU-combined-ali-a2p-RAXML-GTRGAMMA-100xRAPIDBP -m GTRGAMMA -N 100 -f a -x 1975021703574 -p 197502170357 

    #Delete long branch
    cp RAxML_bipartitions.SSU-combined-ali-a2p-RAXML-GTRGAMMA-100xRAPIDBP raxml.tre
    run_treeshrink.py -t raxml.tre
    grep ">" Outall_16s-nodupe.fna | cut -f2 -d'>' >  Out-GCA.list
    cp raxml_treeshrink/output_summary.txt .
    $companion output_summary.txt --mode=shrink --shrinkvalue=$shrink
    cp RAxML_bipartitions.SSU-combined-ali-a2p-RAXML-GTRGAMMA-100xRAPIDBP Constained-tree.tre
    prune-tree.pl --negate-list Constained-tree.idl 

    #format RAxlm output mapping idm file
    #refgroup part
    fetch-tax.pl GC.list --taxdir=$taxdir --item-type=taxid
    cut -f1 GC.tax > f1.tax
    cut -f2 GC.tax > f2.tax
    paste f2.tax f1.tax > part1.idm
    #SSU User file part
    grep ">" SSU.fasta | cut -f2 -d ">" | cut -f1 -d" " > ssu.f1
    grep ">" SSU.fasta | cut -f2 -d ">"  > ssu.f2
    paste ssu.f2 ssu.f1 > part2.idm
    #combined
    cat part1.idm part2.idm > SSU.idm
    cp SSU.idm RAXML.idm
    mv Constained-tree.tre RAXML.tre
    format-tree.pl RAXML.tre --map-ids --figtree
    format-tree.pl RAXML.tre --map-ids
    mv RAXML.nex Constained-tree.nex
    mv RAXML.tre Constained-tree.tre
    """
}

//output the results
process PublicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'Constained-tree.nex' from constTreeNexus_ch
    file 'Constained-tree.tre' from constTreeTre_ch
    file 'SSU-combined-ali-a2p.fasta' from alignementBmge_ch
    file 'SSU-combined-ali.fasta' from alignement_ch
    file 'GC.list' from gcfList_ch

    output:
    file 'Constained-tree.nex' into outConstTreeNexus_ch
    file 'Constained-tree.tre' into outConstTreeTre_ch
    file 'SSU-combined-ali-a2p.fasta' into outAlignementBmge_ch
    file 'SSU-combined-ali.fasta' into outAlignement_ch
    file 'GC.list' into outGcfList_ch

    //script
    script:

    """
    cp Constained-tree.tre Constained-tree-FINAL.tre
    cp Constained-tree.nex Constained-tree-FINAL.nex
    cp SSU-combined-ali-a2p.fasta SSU-combined-ali-a2p-FINAL.fasta
    cp SSU-combined-ali.fasta SSU-combined-ali-FINAL.fasta
    cp GC.list GC-FINAL.list
    echo "ORPER analyses completed" > log
    """
}





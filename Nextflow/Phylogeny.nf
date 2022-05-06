#!/usr/bin/env nextflow

/*
========================================================================================
                         Phylogeny GENERA
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

    Version: 1.0.2 
    
    Citation:
    Please cite : 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow Phylogeny.nf --OG=OGs --IDM=file.idm --jackk=yes 
    
    Mandatory arguments:
    --OG                    Path to OG directory in fasta format (.faa for prot and .fna files for DNA)
    --IDM                   Path to IDM file

    Optional arguments:
    --mode                  specify prot or DNA, default = prot
    --align                 activate the alignment for protein, yes or no, default = no
    --jackk                 activate jackknife, default = no
    --rep                   number of jackkninfe replicates
    --width                 Width of jackknife replicates
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

//Path to OGs : Mandatory
params.OG = null
if (params.OG == null) {
	exit 1, "Path to OG directory in fasta format."
}

//Path to IDM
params.IDM = null
if (params.IDM == null) {
	exit 1, "Path to IDM file containing final names."
}

//mode
params.mode = 'prot'

//align
params.align = 'no'

//jaccknife
params.jackk = 'no'

//Jaccknife replicates number
params.rep = '100'

//Jaccknife width
params.width = '100000'

//outdir
params.outdir='GENERA-phylogeny'

//cpu
params.cpu = '1'

//Path to companion
params.companion = '/opt/companion/Phylogeny_companion.py'

//version
params.version = '1.0.2'

/*
CORE PROGRAM
*/

//Load input files
og_ch = Channel.fromPath(params.OG)
idm_ch = Channel.fromPath(params.IDM)


//change fasta 
process changext {
	//informations

	//input output
    input:
    file '*.f*' from og_ch

    output:
    file 'EXT' into ext_ch1

    //script
    script:
    println "GENERA info: format names with BMC"
    if (params.mode == 'prot') {
        println "GENERA info: fasta to faa"
        """
        mkdir EXT
        cp .f/* EXT
        cd EXT
        find *.* | cut -f1 -d'.' > list
        for f in `cat list`; do mv \$f.fasta \$f.faa; done
        cd ../
        """
    }
    else if (params.mode == 'DNA') {
        println "GENERA info: fasta to fna"
        """
        mkdir EXT
        cp .f/* EXT
        cd EXT
        find *.* | cut -f1 -d'.' > list
        for f in `cat list`; do mv \$f.fasta \$f.fna; done
        cd ../
        """      
    }
}

//BMC names format
process format {
	//informations

	//input output
    input:
    file 'EXT' from ext_ch1

    output:
    file 'FASTA/*abbr*' into abbr_ch1
    file "GENERA-Phylogeny.log" into log_ch1

    //script
    script:
    println "GENERA info: format names with BMC"
    if (params.mode == 'prot') {
        println "GENERA info: Phylogeny protein mode"
        """
        mkdir FASTA
        cp EXT/* FASTA/
        cd FASTA
        find *.* | cut -f1 -d'.' > list
        for f in `cat list`; do sed -i -e 's/@/ /g' \$f*.faa; inst-abbr-ids.pl \$f*.faa --id-regex=:DEF; sed -i -e 's/|//g' \$f*abbr*.faa; done
        cd ../
        echo "GENERA info: Phylogeny protein mode" >> GENERA-Phylogeny.log
        echo "GENERA info: format names with BMC" >> GENERA-Phylogeny.log
        """
    }
    else if (params.mode == 'DNA') {
        println "GENERA info: Phylogeny DNA mode"
        """
        mkdir FASTA
        cp EXT/* FASTA/
        cd FASTA
        find *.* | cut -f1 -d'.' > list
        for f in `cat list`; do sed -i -e 's/@/ /g' \$f*.fna; inst-abbr-ids.pl \$f*.fna --id-regex=:DEF; sed -i -e 's/|//g' \$f*abbr*.fna; done
        cd ../
        echo "GENERA info: Phylogeny DNA mode" >> GENERA-Phylogeny.log
        echo "GENERA info: format names with BMC" >> GENERA-Phylogeny.log
        """      
    }
}


//Alignment
process alignment {
	//informations

	//input output
    input:
    file '*' from abbr_ch1
    file "GENERA-Orthology.log" from log_ch1

    output:
    file 'aligned/*' into ali_ch1
    file "GENERA-Phylogeny.log" into log_ch2

    //script
    script:
    println "GENERA info: run alignment"
    if (params.mode == 'prot') {
        if (params.align == 'yes') {
            """
            mkdir aligned
            mkdir OGs
            find *.faa | cut -f1 -d"." > list
            mv *.faa OGs/
            for f in `cat list`; do muscle3.8.31_i86linux64 -in OGs/\$f.faa -out aligned/\$f.faa; done
            echo "GENERA info: run prot alignment" >> GENERA-Phylogeny.log
            """
        }
        else if (params.align == 'no') {
            """
            mkdir aligned
            mkdir OGs
            find *.faa | cut -f1 -d"." > list
            mv *.faa aligned/
            echo "GENERA info:  No prot alignment" >> GENERA-Phylogeny.log
            """
        }
    }
    else if (params.mode == 'DNA') {
        """
        mkdir aligned
        mv *.fna aligned/
        echo "GENERA info: No alignment for DNA files" >> GENERA-Phylogeny.log
        """
    }
}

//Unambiguous position selection
process unambiguousPosition {
	//informations

	//input output
    input:
    file '*' from ali_ch1
    file "GENERA-Phylogeny.log" from log_ch2

    output:
    file 'UNAMBIG/*.ali' into unambigous_ch1
    file 'UNAMBIG/*.ali' into unambigous_ch2
    file "GENERA-Phylogeny.log" into log_ch3

    //script
    script:
    println "GENERA info: Unambiguous position conservation"
    if (params.mode == 'prot') {
        """
        mkdir BMGE
        mkdir UNAMBIG
        mv *.faa BMGE/
        cd BMGE
        fasta2ali.pl *.faa
        ali2phylip.pl *.ali --bmge-mask=medium --ali
        cd ../
        mv BMGE/*a2p.ali UNAMBIG/
        echo "GENERA info: Unambiguous position conservation" >> GENERA-Phylogeny.log
        """
    }
    else if (params.mode == 'DNA') {
        """
        mkdir UNAMBIG
        mv *.fna UNAMBIG
        cd UNAMBIG/
        fasta2ali.pl *.fna
        cd ../
        echo "GENERA info: No unambiguous position conservation for DNA files" >> GENERA-Phylogeny.log
        """
    }
}

//Concatenation
process concatenation {
	//informations

	//input output
    input:
    file '*' from unambigous_ch1
    val companion from params.companion
    file "GENERA-Phylogeny.log" from log_ch3

    output:
    file 'data-ass/data-ass.ali' into concat_ch1
    file "GENERA-Phylogeny.log" into log_ch4

    //script
    script:
    println "GENERA info: Concatenation"
    if (params.mode == 'prot') {
        """
        echo 'runniner iter' > iter
        $companion iter --mode=iter
        mkdir ali
        mkdir raw/
        mv *-iter.ali ali/
        mv *.ali raw/
        scafos in=ali out=otu
        scafos in=ali out=data otu=otu/otu-freq.otu o=ov
        mkdir final-ali/
        mv data/*.ali final-ali/
        scafos in=final-ali/ out=data-ass otu=otu/otu-freq.otu gamma=yes o=gclv g=30 format=fpm
        echo "GENERA info: Concatenation" >> GENERA-Phylogeny.log
        """
    }
    else if (params.mode == 'DNA') {
        """
        echo 'runniner iter' > iter
        $companion iter --mode=iter
        mkdir ali
        mkdir raw/
        mv *-iter.ali ali/
        mv *.ali raw/
        #mv *.ali ali/
        scafos in=ali out=otu
        scafos in=ali out=data otu=otu/otu-freq.otu o=ov
        mkdir final-ali/
        mv data/*.ali final-ali/
        scafos in=final-ali/ out=data-ass otu=otu/otu-freq.otu gamma=yes o=gclv g=30 format=fpm
        echo "GENERA info: Concatenation" >> GENERA-Phylogeny.log
        """
    }
}


//ML inference
process mlinference {
	//informations

	//input output
    input:
    file 'data-ass.ali' from concat_ch1
    file "GENERA-Phylogeny.log" from log_ch4
    val cpu from params.cpu
    val companion from params.companion

    output:
    file 'protML.tre' into mltree_ch1
    file 'DNAclassic.tre' into mlDNAclassic_ch1
    file 'DNAtwo.tre' into mlDNAtwo_ch1
    file 'DNAthird.tre' into mlDNAthird_ch1
    file "GENERA-Phylogeny.log" into log_ch5

    //script
    script:
    if (params.mode == 'prot') {
        println "GENERA info: Prot ML inference"
        """
        ali2phylip.pl data-ass.ali --map-ids
        raxmlHPC-PTHREADS-AVX -T $cpu -s data-ass.phy \
        -n data-ass-RAXML-PROTGAMMALGF-100xRAPIDBP -m PROTGAMMALGF -N 100 -f a -x 1975021703574 -p 1975021703574
        #format
        mv data-ass.idm RAxML_bipartitions.idm; format-tree.pl RAxML_bipartitions.data-ass-RAXML-PROTGAMMALGF-100xRAPIDBP --map-ids; mv RAxML_bipartitions.tre protML.tre
        echo "GENERA info: ML inference" >> GENERA-Phylogeny.log
        #False file part(for DNA)
        echo 'FALSE IDM' > data-ass-blocks.idm
        echo 'FALSE DNA RAXML' > DNAclassic.tre
        echo 'FALSE DNA RAXML' > DNAtwo.tre
        echo 'FALSE DNA RAXML' > DNAthird.tre
        """
    }
    else if (params.mode == 'DNA') {
        println "GENERA info: DNA ML inference"
        """
        #Matrix phylip convert
        ali2phylip.pl data-ass.ali --map-ids
        #RaxML
        #Without codon partition
        raxmlHPC-PTHREADS-AVX -T $cpu -s data-ass.phy \
        -n classic -m GTRGAMMA -N 100 -f a -x 1975021703574 -p 1975021703574
        cp data-ass.idm RAxML_bipartitions.idm; format-tree.pl RAxML_bipartitions.classic --map-ids; mv RAxML_bipartitions.tre DNAclassic.tre
        echo "GENERA info: DNA ML inference, no partition" >> GENERA-Phylogeny.log

        #Two first codon only
        $companion data-ass.phy --mode=two
        mask-ali.pl blocks --alifile=data-ass.ali
        ali2phylip.pl data-ass-blocks.ali --map-ids
        raxmlHPC-PTHREADS-AVX -T $cpu -s data-ass-blocks.phy \
        -n two -m GTRGAMMA -N 100 -f a -x 1975021703574 -p 1975021703574
        cp data-ass-blocks.idm RAxML_bipartitions.idm; format-tree.pl RAxML_bipartitions.two --map-ids; mv RAxML_bipartitions.tre DNAtwo.tre
        echo "GENERA info: DNA ML inference, codon 1&2" >> GENERA-Phylogeny.log           

        #Partition on third codon
        $companion data-ass.phy --mode=third
        raxmlHPC-PTHREADS-AVX -T $cpu -s data-ass.phy -m GTRGAMMA -q partition.txt \
        -n third -N 100 -f a -x 1975021703574 -p 1975021703574
        cp data-ass.idm RAxML_bipartitions.idm; format-tree.pl RAxML_bipartitions.third --map-ids; mv RAxML_bipartitions.tre DNAthird.tre
        echo "GENERA info: DNA ML inference, third codon partition" >> GENERA-Phylogeny.log

        #False file part(for Prot)
        echo 'FALSE PROT RAXML' > protML.tre
        """
    }
}

//Jackknife matrices
process jaccknifeMatrix {
	//informations

	//input output
    input:
    file '*' from unambigous_ch2
    file "GENERA-Phylogeny.log" from log_ch5
    val rep from params.rep
    val width from params.width

    output:
    file 'matrices/*.ali' into jackkconcats_ch1
    file "GENERA-Phylogeny.log" into log_ch6

    //script
    script:
    if (params.mode == 'prot') {
        if (params.jackk == 'yes') {
            println "GENERA info: make jaccknife matrices"
            """
            #construct datasets
            mkdir ali
            mv *.ali ali/
            jack-ali-dir.pl ali/ --rep=$rep --width=$width
            #Go into each replicate and generate supermatrices
            cd ali-jack*/
            mv replicate-* ../
            cd ../
            ls -d replicate*/ | cut -f1 -d'/' > list
            for f in `cat list`; do cd \$f; mkdir ali; mv *.ali ali/; scafos in=ali out=otu; scafos in=ali out=data otu=otu/otu-freq.otu o=ov; \
            mkdir final-ali/; mv data/*.ali final-ali/; scafos in=final-ali/ out=data-ass otu=otu/otu-freq.otu gamma=yes o=gclv g=30 format=fpm; cd  ../; done
            #Storage of matrices
            mkdir matrices
            for f in `cat list`; do cd \$f; cd data-ass/; mv data-ass.ali ../../matrices/data-ass_\$f.ali; cd ../../; done
            echo "GENERA info: make jaccknife matrices" >> GENERA-Phylogeny.log
            """
        }
        else {
            println "GENERA info: jaccknife not activated"
            """
            mkdir matrices
            echo "GENERA info: jaccknife not activated" > matrices/FALSE-jackk.ali
            echo "GENERA info: jaccknife not activated" >> GENERA-Phylogeny.log
            """
        }
    }
    else if (params.mode == 'DNA') {
        if (params.jackk == 'yes') {
            println "GENERA info: make jaccknife matrices"
            """
            #construct datasets
            mkdir ali
            mv *.ali ali/
            jack-ali-dir.pl ali/ --rep=$rep --width=$width
            #Go into each replicate and generate supermatrices
            cd ali-jack*/
            mv replicate-* ../
            cd ../
            ls -d replicate*/ | cut -f1 -d'/' > list
            for f in `cat list`; do cd \$f; mkdir ali; mv *.ali ali/; scafos in=ali out=otu; scafos in=ali out=data otu=otu/otu-freq.otu o=ov; \
            mkdir final-ali/; mv data/*.ali final-ali/; scafos in=final-ali/ out=data-ass otu=otu/otu-freq.otu gamma=yes o=gclv g=30 format=fpm; cd  ../; done
            #Storage of matrices
            mkdir matrices
            for f in `cat list`; do cd \$f; cd data-ass/; mv data-ass.ali ../../matrices/data-ass_\$f.ali; cd ../../; done
            echo "GENERA info: make jaccknife matrices" >> GENERA-Phylogeny.log
            """            
        }
        else {
            println "GENERA info: jaccknife not activated"
            """
            mkdir matrices
            echo "GENERA info: jaccknife not activated" > matrices/FALSE-jackk.ali
            echo "GENERA info: jaccknife not activated" >> GENERA-Phylogeny.log
            """           
        }
    }
}

//Jackknife ML inference
process jaccknifeMLinference {
	//informations

	//input output
    input:
    file '*' from jackkconcats_ch1
    file "GENERA-Phylogeny.log" from log_ch6
    val cpu from params.cpu
    val companion from params.companion

    output:
    file 'trees/jackk.tre' into jackktre_ch1
    file 'trees/jackk-two.tre' into jackktretwo_ch1
    file 'trees/jackk-third.tre' into jackktrethird_ch1
    file "GENERA-Phylogeny.log" into log_ch7

    //script
    script:
    if (params.mode == 'prot') {
        if (params.jackk == 'yes') {
            println "GENERA info: jaccknife Prot ML inferences"
            """
            #format matrices
            find *.ali | cut -f1 -d'.' > list
            for f in `cat list`; do mkdir \$f; mv \$f.ali \$f/; done
            for f in `cat list`; do cd \$f; ali2phylip.pl \$f.ali --map-ids; cd ../; done
            #ML inference
            for f in `cat list`; do cd \$f; raxmlHPC-PTHREADS-AVX -T $cpu -s \$f.phy -n \$f-RAXML-PROTGAMMALGF-fast -m PROTGAMMALGF -p 1975021703574  -f F; cd ../; done
            #format tree
            mkdir protJ
            for f in `cat list`; do cd \$f; mv *.idm temp.idm; cut -f1 temp.idm > f1; cut -f2 temp.idm > f2; sed -i -e 's/ /_/g' f1; sed -i -e 's/#//g' f2; cut -f1 -d"@" f1 > temp; \
            mv temp f1; paste f1 f2 > RAxML_fastTree.idm;  format-tree.pl RAxML_fastTree.*-RAXML-PROTGAMMALGF-fast --map-ids; \
            mv RAxML_fastTree.tre ../protJ/\$f.tre; cd ../; done   
            #make consensus tree
            cd protJ
            cat *.tre > intree
            echo "consense << EOF" >> makecons
            echo '\$1' >> makecons
            echo "y" >> makecons
            echo "EOF" >> makecons
            chmod a+x makecons
            ./makecons intree
            mv outtree jackk.tre
            cd  ../
            mkdir trees
            mv protJ/jackk.tre trees/
            echo "GENERA info: jaccknife ML inferences" >> GENERA-Phylogeny.log
            #False tre for DNA
            echo "GENERA info: jaccknife not activated" > trees/jackk-two.tre
            echo "GENERA info: jaccknife not activated" > trees/jackk-third.tre           
            """
        }
        else {
            println "GENERA info: jaccknife not activated"
            """
            mkdir trees
            echo "GENERA info: jaccknife not activated" > trees/jackk.tre
            echo "GENERA info: jaccknife not activated" > trees/jackk-two.tre
            echo "GENERA info: jaccknife not activated" > trees/jackk-third.tre
            echo "GENERA info: jaccknife not activated" >> GENERA-Phylogeny.log
            """
        }
    }
    else if (params.mode == 'DNA') {
        if (params.jackk == 'yes') {
            println "GENERA info: jaccknife DNA ML inferences"
            """
            #Jackknife without partition
            #format matrices
            find *.ali | cut -f1 -d'.' > list
            for f in `cat list`; do mkdir \$f-class; cp \$f.ali \$f-class/; done
            for f in `cat list`; do cd \$f-class; ali2phylip.pl \$f.ali --map-ids; cd ../; done
            #ML inference
            for f in `cat list`; do cd \$f-class; raxmlHPC-PTHREADS-AVX -T $cpu -s \$f.phy -n \$f-classic-fast -m GTRGAMMA -p 1975021703574  -f F; cd ../; done
            #format tree
            mkdir dnaJclassic
            for f in `cat list`; do cd \$f-class; mv *.idm temp.idm; cut -f1 temp.idm > f1; cut -f2 temp.idm > f2; sed -i -e 's/ /_/g' f1; sed -i -e 's/#//g' f2; \
            cut -f1 -d"@" f1 > temp; mv temp f1; paste f1 f2 > RAxML_fastTree.idm; format-tree.pl RAxML_fastTree.*-classic-fast --map-ids; \
            mv RAxML_fastTree.tre ../dnaJclassic/\$f.tre; cd ../; done   
            #make consensus tree
            cd dnaJclassic
            cat *.tre > intree
            echo "consense << EOF" >> makecons
            echo '\$1' >> makecons
            echo "y" >> makecons
            echo "EOF" >> makecons
            chmod a+x makecons
            ./makecons intree
            mv outtree jackk.tre
            cd  ../
            mkdir trees
            mv dnaJclassic/jackk.tre trees/
            echo "GENERA info: jaccknife DNA inferences, no partition" >> GENERA-Phylogeny.log

            #Jackknife on two first codons
            #format matrices
            find *.ali | cut -f1 -d'.' > list
            for f in `cat list`; do mkdir \$f-two; cp \$f.ali \$f-two/; done
            for f in `cat list`; do cd \$f-two; ali2phylip.pl \$f.ali --map-ids; $companion \$f.phy --mode=two; mask-ali.pl blocks --alifile=\$f.ali; ali2phylip.pl \$f-blocks.ali --map-ids; cd ../; done
            #ML inference
            for f in `cat list`; do cd \$f-two; raxmlHPC-PTHREADS-AVX -T $cpu -s \$f-blocks.phy -n \$f-two-fast -m GTRGAMMA -p 1975021703574  -f F; cd ../; done
            #format tree
            mkdir dnaJtwo
            for f in `cat list`; do cd \$f-two; mv *blocks.idm temp.idm; cut -f1 temp.idm > f1; cut -f2 temp.idm > f2; sed -i -e 's/ /_/g' f1; sed -i -e 's/#//g' f2; \
            cut -f1 -d"@" f1 > temp; mv temp f1; paste f1 f2 > RAxML_fastTree.idm; format-tree.pl RAxML_fastTree.*-two-fast --map-ids; \
            mv RAxML_fastTree.tre ../dnaJtwo/\$f.tre; cd ../; done   
            #make consensus tree
            cd dnaJtwo
            cat *.tre > intree
            echo "consense << EOF" >> makecons
            echo '\$1' >> makecons
            echo "y" >> makecons
            echo "EOF" >> makecons
            chmod a+x makecons
            ./makecons intree
            mv outtree jackk-two.tre
            cd  ../
            #mkdir trees
            mv dnaJtwo/jackk-two.tre trees/
            echo "GENERA info: jaccknife ML inferences, codon 1&2" >> GENERA-Phylogeny.log  

            #Jackknife with partition on third codon
            #format matrices
            find *.ali | cut -f1 -d'.' > list
            for f in `cat list`; do mkdir \$f-third; mv \$f.ali \$f-third/; done
            for f in `cat list`; do cd \$f-third; ali2phylip.pl \$f.ali --map-ids; cd ../; done
            #ML inference
            for f in `cat list`; do cd \$f-third; $companion \$f.phy --mode=third; raxmlHPC-PTHREADS-AVX -T $cpu -s \$f.phy -n \$f-third-fast -m GTRGAMMA -q partition.txt -p 1975021703574  -f F; cd ../; done
            #format tree
            mkdir dnaJthird
            for f in `cat list`; do cd \$f-third; mv *.idm temp.idm; cut -f1 temp.idm > f1; cut -f2 temp.idm > f2; sed -i -e 's/ /_/g' f1; sed -i -e 's/#//g' f2; \
            cut -f1 -d"@" f1 > temp; mv temp f1; paste f1 f2 > RAxML_fastTree.idm; format-tree.pl RAxML_fastTree.*-third-fast --map-ids; \
            mv RAxML_fastTree.tre ../dnaJthird/\$f.tre; cd ../; done   
            #make consensus tree
            cd dnaJthird
            cat *.tre > intree
            echo "consense << EOF" >> makecons
            echo '\$1' >> makecons
            echo "y" >> makecons
            echo "EOF" >> makecons
            chmod a+x makecons
            ./makecons intree
            mv outtree jackk-third.tre
            cd  ../
            #mkdir trees
            mv dnaJthird/jackk-third.tre trees/
            echo "GENERA info: jaccknife DNA inferences, partition on third codon" >> GENERA-Phylogeny.log
            """
        }
        else {
            println "GENERA info: jaccknife not activated"
            """
            mkdir trees
            echo "GENERA info: jaccknife not activated" > trees/jackk.tre
            echo "GENERA info: jaccknife not activated" > trees/jackk-two.tre
            echo "GENERA info: jaccknife not activated" > trees/jackk-third.tre
            echo "GENERA info: jaccknife not activated" >> GENERA-Phylogeny.log
            """
        }
    }
}

//format trees
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'protML.tre' from mltree_ch1
    file 'DNAclassic.tre' from mlDNAclassic_ch1
    file 'DNAtwo.tre' from mlDNAtwo_ch1
    file 'DNAthird.tre' from mlDNAthird_ch1
    file 'jackk.tre' from jackktre_ch1
    file 'jackk-two.tre' from jackktretwo_ch1
    file 'jackk-third.tre' from jackktrethird_ch1
    file 'IDM' from idm_ch
    file "GENERA-Phylogeny.log" from log_ch7
    val companion from params.companion

    output:
    file 'protML-final.tre' into protMLFINAL_ch
    file 'jackk-Prot-final.tre' into jackkprotFINAL_ch
    file 'DNAclassic-final.tre' into dnaclassicFINAL_ch
    file 'DNAtwo-final.tre' into dnatwoFINAL_ch
    file 'DNAthird-final.tre' into dnathirdFINAL_ch
    file 'jackk-DNA-final.tre' into jackkDNAFINAL_ch
    file 'jackk-DNAtwo-final.tre' into jackkDNAtwoFINAL_ch
    file 'jackk-DNAthird-final.tre' into jackkDNAthirdFINAL_ch
    file "GENERA-Phylogeny.log" into logFINAL_ch

    //script
    script:
    println "GENERA info: format-trees"
    if (params.mode == 'prot') {
        if (params.jackk == 'yes') {
            """
            #Prot
            #Large ML
            tree2list.pl protML.tre
            sed -i -e 's/ /_/g' protML.idl
            $companion protML.idl --mode=IDM
            cp idm.temp protML.idm; format-tree.pl protML.tre --map-ids; mv protML.tre protML-final.tre
            #Jack
            cp IDM jackk.idm; format-tree.pl jackk.tre --map-ids --from-consense; mv jackk.tre jackk-Prot-final.tre
            #False DNA
            echo "GENERA info: false file" > DNAclassic-final.tre
            echo "GENERA info: false file" > DNAtwo-final.tre
            echo "GENERA info: false file" > DNAthird-final.tre
            echo "GENERA info: false file" > jackk-DNA-final.tre
            echo "GENERA info: false file" > jackk-DNAtwo-final.tre
            echo "GENERA info: false file" > jackk-DNAthird-final.tre
            echo "GENERA info: format-trees" >> GENERA-Phylogeny.log
            """
        }
        else {
            """
            #Prot
            #Large ML
            tree2list.pl protML.tre
            sed -i -e 's/ /_/g' protML.idl
            $companion protML.idl --mode=IDM
            cp idm.temp protML.idm; format-tree.pl protML.tre --map-ids; mv protML.tre protML-final.tre
            #Jack
            mv jackk.tre jackk-final.tre           
            #False DNA
            echo "GENERA info: false file" > DNAclassic-final.tre
            echo "GENERA info: false file" > DNAtwo-final.tre
            echo "GENERA info: false file" > DNAthird-final.tre
            echo "GENERA info: false file" > jackk-DNA-final.tre
            echo "GENERA info: false file" > jackk-DNAtwo-final.tre
            echo "GENERA info: false file" > jackk-DNAthird-final.tre
            echo "GENERA info: format-trees" >> GENERA-Phylogeny.log
            echo "GENERA info: false file" > jackk-Prot-final.tre
            """
        }
    }
    else if (params.mode == 'DNA') {
        if (params.jackk == 'yes') {
            """
            #DNA
            #Large
            tree2list.pl DNAclassic.tre
            sed -i -e 's/ /_/g' DNAclassic.idl
            $companion DNAclassic.idl --mode=IDM
            cp idm.temp DNAclassic.idm; format-tree.pl DNAclassic.tre --map-ids; mv DNAclassic.tre DNAclassic-final.tre
            tree2list.pl DNAtwo.tre
            sed -i -e 's/ /_/g' DNAtwo.idl
            $companion DNAtwo.idl --mode=IDM           
            cp idm.temp DNAtwo.idm; format-tree.pl DNAtwo.tre --map-ids; mv DNAtwo.tre DNAtwo-final.tre
            tree2list.pl DNAthird.tre
            sed -i -e 's/ /_/g' DNAthird.idl
            $companion DNAthird.idl --mode=IDM
            cp idm.temp DNAthird.idm; format-tree.pl DNAthird.tre --map-ids; mv DNAthird.tre DNAthird-final.tre
            #Jack
            cp IDM jackk.idm; format-tree.pl jackk.tre --map-ids --from-consense; mv jackk.tre jackk-DNA-final.tre  
            cp IDM jackk-two.idm; format-tree.pl jackk-two.tre --map-ids --from-consense; mv jackk-two.tre jackk-DNAtwo-final.tre
            cp IDM jackk-third.idm; format-tree.pl jackk-third.tre --map-ids --from-consense; mv jackk-third.tre jackk-DNAthird-final.tre
            #False prot
            echo "GENERA info: false file" > protML-final.tre
            echo "GENERA info: false file" > jackk-Prot-final.tre        
            echo "GENERA info: format-trees" >> GENERA-Phylogeny.log
            """
        }
        else {
            """
            #Large
            tree2list.pl DNAclassic.tre
            sed -i -e 's/ /_/g' DNAclassic.idl
            $companion DNAclassic.idl --mode=IDM
            cp idm.temp DNAclassic.idm; format-tree.pl DNAclassic.tre --map-ids; mv DNAclassic.tre DNAclassic-final.tre
            tree2list.pl DNAtwo.tre
            sed -i -e 's/ /_/g' DNAtwo.idl
            $companion DNAtwo.idl --mode=IDM           
            cp idm.temp DNAtwo.idm; format-tree.pl DNAtwo.tre --map-ids; mv DNAtwo.tre DNAtwo-final.tre
            tree2list.pl DNAthird.tre
            sed -i -e 's/ /_/g' DNAthird.idl
            $companion DNAthird.idl --mode=IDM
            cp idm.temp DNAthird.idm; format-tree.pl DNAthird.tre --map-ids; mv DNAthird.tre DNAthird-final.tre
            #Jack
            mv jackk.tre jackk-DNA-final.tre  
            mv jackk-two.tre jackk-DNAtwo-final.tre
            mv jackk-third.tre jackk-DNAthird-final.tre
            #False prot
            echo "GENERA info: false file" > protML-final.tre
            echo "GENERA info: false file" > jackk-Prot-final.tre 
            echo "GENERA info: format-trees" >> GENERA-Phylogeny.log
            """
        }
    }

}
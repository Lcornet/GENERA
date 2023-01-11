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

    Version: 2.0.0 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow run ANI.nf --genome=genome --list=list --cpu=20
    
    Mandatory arguments:
    --genome                 Specify the path to genome directory (ext = .fna)
    --list                   Specify the order of organism for the heatmap.
    --idm                    Specify the IDM file

    Optional arguments:
    --tool                   Specify which ANI tool to use: fastANI, orthoani, default == fastANI
    --mode                   Specifiy the mode, onetomany or manytomany, default == manytomany
    --tree                   Activate the tree inference, yes or no, default == no
    --shortlist              Specify the list of query genomes for onetomany mode
    --minFraction            Specify the minimal fraction between two genomes to consider ANI for fastANI, default == 0.2
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
	exit 1, "Specify the path to genome directory (ext = .fna)"
}

//Path to OGs : Mandatory
params.list = null
if (params.list == null) {
	exit 1, "Specify the order of organism for the heatmap."
}

//Path to idm idl
params.idm = null
if (params.idm == null) {
	exit 1, "Specify the order of organism for the heatmap."
}

//tool
params.tool = 'fastANI'

//mode
params.mode = 'manytomany'

//tree
params.tree = 'no'

//short list
params.shortlist = 'none'

//minFraction
params.minFraction = '0.2'

//cpu
params.cpu = '1'

//outdir
params.outdir='GENERA_ANI'

//Path to companion
params.companion = '/opt/companion/ANI_companion.py'

//version
params.version = '2.0.0'

/*
CORE PROGRAM
*/

//Load input files
genome_ch1 = Channel.fromPath(params.genome)
genome_ch2 = Channel.fromPath(params.genome)
list_ch = Channel.fromPath(params.list)
idm_ch1 = Channel.fromPath(params.idm)
idm_ch2 = Channel.fromPath(params.idm)
shortlist_ch1 = Channel.fromPath(params.shortlist)
shortlist_ch2 = Channel.fromPath(params.shortlist)
shortlist_ch3 = Channel.fromPath(params.shortlist)

//BMC names format
process fastANI {
	//informations

	//input output
    input:
    file '*' from genome_ch1
    file 'shortlist' from shortlist_ch1
    val cpu from params.cpu
    val fraction from params.minFraction

    output:
    file 'fANI.txt' into ffile_ch1
    file 'fANI.txt' into ffile_ch2
    file 'fANI.txt' into ffile_ch3

    //script
    script:
    if (params.tool == 'fastANI'){
        if (params.mode == 'onetomany') {
            println "GENERA info: running fastANI"
            """
            mkdir GEN
            cp genome/* GEN/
            cd GEN/
            find *.fna > list
            cp ../shortlist .
            for f in `cat shortlist`; do fastANI -q \$f.fna --rl list -o ANI-\$f -t $cpu --matrix --minFraction fraction ; done
            rm -f *.matrix
            cd ../
            cat GEN/ANI* > fANI.txt 
            """
        }
        else if (params.mode == 'manytomany') {
            println "GENERA info: running fastANI"
            """
            mkdir GEN
            cp genome/* GEN/
            cd GEN/
            find *.fna > list
            fastANI --ql list --rl list -o ANI -t $cpu --matrix --minFraction fraction
            #fastANI --ql list  --rl list -o fastANI.output.txt
            cd ../
            cp GEN/ANI fANI.txt
            """
        }
    }
    else {
        println "GENERA info: NOT running fastANI"
        """
        echo "NOT running fastANI" > fANI.txt
        """
    }
}

//BMC names format
process orthoani {
	//informations

	//input output
    input:
    file '*' from genome_ch2
    file 'shortlist' from shortlist_ch2
    val cpu from params.cpu

    output:
    file 'oANI.txt' into ofile_ch1
    file 'oANI.txt' into ofile_ch2
    file 'oANI.txt' into ofile_ch3


    //script
    script:
    if (params.tool == 'orthoani'){
        if (params.mode == 'onetomany') {
            println "GENERA info: running fastANI"
            """
            mkdir GEN
            cp genome/* GEN/
            cd GEN/
            find *.fna > list
            cp ../shortlist .
            sed -i -e 's/.fna//g' list
            sed -i -e 's/.fna//g' shortlist
            for f in `cat shortlist`; do for d in `cat list`; do orthoani -q \$f.fna -r \$d.fna > temp; echo -n \$f-\$d- > \$f-\$d.tempANI; echo -n \$(<temp) >> \$f-\$d.tempANI; echo -X-X >> \$f-\$d.tempANI; sed -i -e 's/-/\t/g' \$f-\$d.tempANI; rm -f temp; done; done
            cat *.tempANI > ANI
            cd ../
            cat GEN/ANI* > oANI.txt 
            """
        }
        else if (params.mode == 'manytomany') {
            println "GENERA info: running fastANI"
            """
            mkdir GEN
            cp genome/* GEN/
            cd GEN/
            find *.fna > list
            sed -i -e 's/.fna//g' list
            for f in `cat list`; do for d in `cat list`; do orthoani -q \$f.fna -r \$d.fna > temp; echo -n \$f-\$d- > \$f-\$d.tempANI; echo -n \$(<temp) >> \$f-\$d.tempANI; echo -X-X >> \$f-\$d.tempANI; sed -i -e 's/-/\t/g' \$f-\$d.tempANI; rm -f temp; done; done
            cat *.tempANI > ANI
            cd ../
            cp GEN/ANI oANI.txt
            """
        }
    }
    else {
        println "GENERA info: NOT running ortho"
        """
        echo "NOT running orthoani" > oANI.txt
        """
    }
}

//modellling plots
process heatmap {
	//informations

	//input output
    input:
    file 'fANI.txt' from ffile_ch1
    file 'oANI.txt' from ofile_ch1
    file 'list' from list_ch
    file 'shortlist' from shortlist_ch3
    file 'file.idm' from idm_ch1
    val companion from params.companion

    output:
    file 'ANI.pdf' into pdf_ch
    file 'ANI-close' into close_ch
    file 'ANI-matrix.txt' into matrice_ch
    file 'ANI-dist-matrix.txt' into distmatrix_ch

    //script
    script:
    if (params.tool == 'fastANI'){
        if (params.list == 'none'){
            """
            echo "GENERA info: a list of file should be provided" >> GENERA-METABO.log
            """
        }
        else if (params.mode == 'onetomany') {
            println "GENERA info: heatmap onetomany mode"
            """
            mv fANI.txt ANI.txt
            $companion ANI.txt --mode=onetomany
            #Don't run R in onetomany mode, jusr closest list
            echo "onetomany mode" > ANI.pdf
            #close part
            mkdir ANI-close
            mv *closest* ANI-close/
            """        
        }
        else {
            println "GENERA info: heatmap, manytomany mode"
            """
            mv fANI.txt ANI.txt
            $companion ANI.txt --mode=manytomany --anitool=fastANI
            $companion ANI.txt --mode=idm --submode=heatmap 
            rm -f ANI-infile.txt
            mv ANI-heat.txt ANI-infile.txt
            #run R
            #cp list finallist
            for f in `cat list`; do grep \$f file.idm | cut -f1 ; done > finallist
            echo "GENOME" > head
            cat head finallist > temp; mv -f temp finallist
            tpage --define file=ANI /opt/code-ANI.tt > code.r
            Rscript code.r
            mv code.r ANI-code.r
            #close part
            mkdir ANI-close
            mv *closest* ANI-close/
            """
        }
    }
    else if (params.tool == 'orthoani') {
        if (params.list == 'none'){
            """
            echo "GENERA info: a list of file should be provided" >> GENERA-METABO.log
            """
        }
        else if (params.mode == 'onetomany') {
            println "GENERA info: heatmap onetomany mode"
            """
            mv oANI.txt ANI.txt
            $companion ANI.txt --mode=onetomany
            #Don't run R in onetomany mode, jusr closest list
            echo "onetomany mode" > ANI.pdf
            #close part
            mkdir ANI-close
            mv *closest* ANI-close/
            """        
        }
        else {
            println "GENERA info: heatmap, manytomany mode"
            """
            mv oANI.txt ANI.txt
            $companion ANI.txt --mode=manytomany --anitool=orthoani
            $companion ANI.txt --mode=idm --submode=heatmap 
            rm -f ANI-infile.txt
            mv ANI-heat.txt ANI-infile.txt
            #run R
            #cp list finallist
            for f in `cat list`; do grep \$f file.idm | cut -f1 ; done > finallist
            echo "GENOME" > head
            cat head finallist > temp; mv -f temp finallist
            tpage --define file=ANI /opt/code-ANI.tt > code.r
            Rscript code.r
            mv code.r ANI-code.r
            #close part
            mkdir ANI-close
            mv *closest* ANI-close/
            """
        }
    }
}

//modellling plots
process NJtree {
	//informations

	//input output
    input:
    file 'fANI.txt' from ffile_ch2
    file 'oANI.txt' from ofile_ch2
    file 'ANI-dist-matrix.txt' from distmatrix_ch
    file 'file.idm' from idm_ch2
    val companion from params.companion

    output:
    file 'ANI.phylip' into distance_ch1
    file 'ANI.newick' into tree_ch1

    //script
    script:
    if (params.mode == 'onetomany') {
        println "GENERA info: heatmap onetomany mode"
        """
        #Onetomany mode, no compute of NJ tree
        echo "Onetomany mode, no compute of NJ tree" > ANI.phylip
        echo "Onetomany mode, no compute of NJ tree" > ANI.newick
        """        
    }
    else {
        if (params.tree == 'yes') {
            println "GENERA info: heatmap, manytomany mode, tree activated"
            """
            #manytomano, inference of NJ tree based on distance matrix
            #mv fANI.txt ANI.txt
            #pairwise_identities_to_distance_matrix.py --max_dist 0.2 ANI.txt > ANI.phylip
            #bionj_tree.R ANI.phylip ANI.newick
            $companion ANI-dist-matrix.txt --mode=idm --submode=dist
            bionj_tree.R ANI.phylip ANI.newick
            """
        }
        else {
            println "GENERA info: heatmap, manytomany mode, tree not activated"
            """
            #Onetomany mode, no compute of NJ tree
            echo "no compute of NJ tree" > ANI.phylip
            echo "no compute of NJ tree" > ANI.newick
            """ 
        }
    }
}


//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'fANI.txt' from ffile_ch3
    file 'oANI.txt' from ofile_ch3
    file 'ANI.pdf' from pdf_ch
    file 'ANI-close' from close_ch
    file 'ANI-matrix.txt' from matrice_ch
    file 'ANI.phylip' from distance_ch1
    file 'ANI.newick' from tree_ch1
    val version from params.version

    output:
    file 'fastANI.txt' into fanifinal_ch
    file 'orthoani.txt' into oanifinal_ch
    file 'ANI.matrix' into matrixfinal_ch
    file 'heatmap.pdf' into heatmapfinal_ch
    file 'ANI-close-match' into closefinal_ch
    file 'ANI-dist.phylip' into distancefinal_ch
    file 'ANI-dist.newick' into treefinal_ch
    file 'GENERA-ANI.log' into logfinal_ch

    //script
    script:
    """
    cp ANI-matrix.txt ANI.matrix
    cp fANI.txt fastANI.txt
    cp oANI.txt orthoani.txt
    cp ANI.pdf heatmap.pdf
    cp -r ANI-close/ ANI-close-match
    cp ANI.phylip ANI-dist.phylip
    cp ANI.newick ANI-dist.newick
    echo "GENERA info: Running fastANI" >> GENERA-ANI.log
    echo VERSION: >> GENERA-ANI.log
    echo $version >> GENERA-ANI.log
    """
}
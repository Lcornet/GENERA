#!/usr/bin/env nextflow

/*
========================================================================================
                         Assembly GENERA
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
    This tool allows to assembly short reads (Illumina) or long reads (Pb/ONT).
    Assembly.nf works for metagenomes or classical assembly.
    Version: 1.0.0 
    
    Citation:
    Please cite : 

    Usage:
    
    The typical command for running the pipeline is as follows:

    nextflow Assembly.nf
    
    Mandatory arguments:
    --shortreadsR1            Path to Illumina short reads forward, R1 fastq.
    --shortreadsR2            Path to Illumina short reads reverse, R2 fastq.                

    Optional arguments:
    --ontreads                Path to Nanopore long reads, fastq
    --pbreads                 Path to Pacbio long reads, fastq
    --genomeSIZE              Size of the expected genome, mandatory for long reads assembler. (m or gb)
    --longreadassembler       Specify which tool to use, CANU or Flye, default = Flye
    --metagenome              metagenome assembly, not activated by default - yes or no.
                              the assembly will be done by metaSPAdes if only illumina short reads are provided
                              if long reads are provided, metaFlye will be used.
    --binner                  Specify which tool to use for binning, metabat or concoct or all,  default = none
    --ragtag                  activate ragat scaffold to map your assembly according to a reference genome (to provide), yes or no, default = no
    --refgenome               reference genome used by ragtag, mandatory with ragtag, default = none
    --outdir                  Name of the output directory, default GENERA-assembly
    --cpu                     Number of cpus, default = 1
    --RAM                     Ram used, default = 500 (Gb)

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

//Path to short reads fastq : Mandatory
params.shortreadsR1 = null
if (params.shortreadsR1 == null) {
	exit 1, "Path to Illumina short reads forward, R1 fastq."
}

//Path to short reads fastq : Mandatory
params.shortreadsR2 = null
if (params.shortreadsR2 == null) {
	exit 1, "Path to Illumina short reads forward, R2 fastq."
}

//Path to nanopore reads fastq
params.ontreads = 'no'

//Path to nanopore reads fastq
params.pbreads = 'no'

//Path to nanopore reads fastq
params.metagenome = 'no'

//Number of cpus
params.cpu = '1'

//RAM used by spades
params.RAM = '280'

//Genome size: Mandatory if Long reads
params.genomeSIZE = null

//Long read assembler
params.longreadassembler= 'Flye'

//binning tool
params.binner = 'none'

//ragoo
params.ragtag = 'no'

//Reference genome used by Ragoo
params.refgenome = 'no'

//outdir
params.outdir='GENERA-assembly'

//Path to companion
params.companion = '/opt/companion/Assembly_companion.py'

//version
params.version = '1.0.0'

/*
CORE PROGRAM
*/

//Load input files
shortreadsR1_ch = Channel.fromPath(params.shortreadsR1)
shortreadsR2_ch = Channel.fromPath(params.shortreadsR2)
ontreads_ch1 = Channel.fromPath(params.ontreads)
ontreads_ch2 = Channel.fromPath(params.ontreads)
pbreads_ch1 = Channel.fromPath(params.pbreads)
pbreads_ch2 = Channel.fromPath(params.pbreads)
refgenome_ch = Channel.fromPath(params.refgenome)

//Trimming of short reads
//Download RefSeq metadata, compute taxonomy file and prudce download ftp file
process Trimming {
	//informations

	//input output
    input:
    file 'R1.fastq' from shortreadsR1_ch
    file 'R2.fastq' from shortreadsR2_ch
    val cpu from params.cpu

    output:
    file "R1-fastp.fastq" into fastpR1_ch1
    file "R2-fastp.fastq" into fastpR2_ch1
    file "R1-fastp.fastq" into fastpR1_ch2
    file "R2-fastp.fastq" into fastpR2_ch2  
    file "R1-fastp.fastq" into fastpR1_ch3
    file "R2-fastp.fastq" into fastpR2_ch3 
    file "R1-fastp_fastqc.html" into fastqcR1_ch
    file "R2-fastp_fastqc.html" into fastqcR2_ch
    file "fastp.log" into logfastp_ch
    file "GENERA-Assembler.log" into log_ch1

    //script
    script:

    """
    fastp --in1 R1.fastq --in2 R2.fastq --out1 R1-fastp.fastq --out2 R2-fastp.fastq --thread $cpu
    fastqc R1-fastp.fastq -t $cpu
    fastqc R2-fastp.fastq -t $cpu
    echo "fastp done" > fastp.log
    echo "GENERA info: Trimming done" >> GENERA-Assembler.log
    """
}

//Short reads assembly
process shortreadsassembly {
    //informations

	//input output
    input:
    file "R1-fastp.fastq" from fastpR1_ch1
    file "R2-fastp.fastq" from fastpR2_ch1
    val cpu from params.cpu
    val ram from params.RAM
    file "GENERA-Assembler.log" from log_ch1

    output:
    file 'SPADES_scaffolds.fasta' into spades_ch1
    file "GENERA-Assembler.log" into log_ch2
    
    //script
    script:
    if (params.ontreads != 'no'){
        println "GENERA info: ONT reads provided, skipping SPAdes"
        """
        echo "GENERA info: ONT reads provided, skipping SPAdes" > SPADES_scaffolds.fasta
        echo "GENERA info: ONT reads provided, skipping SPAdes" >> GENERA-Assembler.log
        """
    }
    else if (params.pbreads != 'no'){
    println "GENERA info: Pacbio reads provided, skipping SPAdes"
        """
        echo "GENERA info: Pacbio reads provided, skipping SPAdes" > SPADES_scaffolds.fasta
        echo "GENERA info: Pacbio reads provided, skipping SPAdes" >> GENERA-Assembler.log
        """
    }
    else if (params.metagenome == 'yes'){
    println "GENERA info: metagenome for short reads activated, Assembly will be done with SPAdes"
        """
        spades.py --pe1-1 `ls R1-fastp.fastq` --pe1-2 `ls R1-fastp.fastq` --meta -t $cpu -m $ram -o SPADES
        cp SPADES/scaffolds.fasta SPADES_scaffolds.fasta
        echo "GENERA info: metagenome for short reads activated, Assembly will be done with SPAdes" >> GENERA-Assembler.log
        """
    }
    else{
    println "GENERA info: Assembly will be done with SPAdes"
        """
        spades.py --pe1-1 `ls R1-fastp.fastq` --pe1-2 `ls R1-fastp.fastq` -t $cpu -m $ram -o SPADES
        cp SPADES/scaffolds.fasta SPADES_scaffolds.fasta
        echo "GENERA info: Assembly will be done with SPAdes" >> GENERA-Assembler.log
        """  
    }  
}

//Long reads assembly
process longreadsassembly {
    //informations

	//input output
    input:
    file 'ont.fastq' from ontreads_ch1
    file 'pb.fastq' from pbreads_ch1
    val cpu from params.cpu
    val size from params.genomeSIZE
    val ram from params.RAM
    file "fastp.log" from logfastp_ch
    file "GENERA-Assembler.log" from log_ch2
    
    output:
    file "LR_assembly.fasta" into lr_ch1
    file "GENERA-Assembler.log" into log_ch3
    
    //script
    script:
    if (params.ontreads != 'no'){
        println "GENERA info: ONT reads provided for assembly"
        if (params.metagenome == 'yes'){
            println "GENERA info: metagenome for ONT long reads activated, Assembly will be done by Flye"
            """
            flye --nano-raw ont.fastq --out-dir FLYE --threads $cpu --meta
            cp FLYE/assembly.fasta LR_assembly.fasta
            echo "GENERA info: metagenome for ONT long reads activated, Assembly will be done by Flye" >> GENERA-Assembler.log
            """
        }
        else if (params.longreadassembler == 'Flye'){
            if (params.genomeSIZE == null){
                println "GENERA info: For long reads assembly, genome size is required"
                """
                """
            }
            else {
                println "GENERA info: Assembly for ONT long reads activated, Assembly will be done by Flye"
                """
                flye --nano-raw ont.fastq --out-dir FLYE --genome-size $size --threads $cpu
                cp FLYE/assembly.fasta LR_assembly.fasta
                echo "GENERA info: Assembly for ONT long reads activated, Assembly will be done by Flye with: genome-size = $size " >> GENERA-Assembler.log
                """
            }
        }
        else if (params.longreadassembler == 'CANU'){
            if (params.genomeSIZE == null){
                println "GENERA info: For long reads assembly, genome size is required"
                """
                """
            }
            else {
                println "GENERA info: Assembly for ONT long reads activated, Assembly will be done by CANU"
                """
                canu -nanopore-raw ont.fastq genomeSize=$size maxMemory=$ram maxThreads=$cpu stopOnLowCoverage=5 cnsErrorRate=0.25 -p CANU -d CANU useGrid=false
                cp CANU/CANU.contigs.fasta LR_assembly.fasta
                echo "GENERA info: Assembly for ONT long reads activated, Assembly will be done by CANU with genomeSize=$size; stopOnLowCoverage=5; cnsErrorRate=0.25 " >> GENERA-Assembler.log
                """
            }
        }
    }
    else if (params.pbreads != 'no'){
        println "GENERA info: Pacbio reads provided for assembly"
        if (params.metagenome == 'yes'){
            println "GENERA info: metagenome for PB long reads activated, assembly will be done by Flye"
            """
            flye --pacbio-raw pb.fastq --out-dir FLYE --genome-size $size --threads $cpu --meta
            cp FLYE/assembly.fasta LR_assembly.fasta
            echo "GENERA info: metagenome for PB long reads activated, Assembly will be done by Flye" >> GENERA-Assembler.log
            """
        }
        else if (params.longreadassembler == 'Flye'){
            if (params.genomeSIZE == null){
                println "GENERA info: For long reads assembly, genome size is required"
                """
                """
            }
            else {
                println "GENERA info: Assembly for PB long reads activated, Assembly will be done by Flye"
                """
                flye --pacbio-raw pb.fastq --out-dir FLYE --genome-size $size --threads $cpu
                cp FLYE/assembly.fasta LR_assembly.fasta
                echo "GENERA info: Assembly for PB long reads activated, Assembly will be done by Flye with: genome-size=$size " >> GENERA-Assembler.log
                """
            }
        }
        else if (params.longreadassembler == 'CANU'){
            if (params.genomeSIZE == null){
                println "GENERA info: For long reads assembly, genome size is required"
                """
                """
            }
            else {
                """
                canu -pacbio-raw pb.fastq genomeSize=$size maxMemory=$ram maxThreads=$cpu stopOnLowCoverage=5 cnsErrorRate=0.25 -p CANU -d CANU useGrid=false
                cp CANU/CANU.contigs.fasta LR_assembly.fasta
                echo "GENERA info: Assembly for PB long reads activated, Assembly will be done by CANU with genomeSize=$size; stopOnLowCoverage=5; cnsErrorRate=0.25 " >> GENERA-Assembler.log
                """
            }
        }
    }
    else {
        println "GENERA info: No long reads provided, skipping long reads assembly"
        """
        echo "GENERA info: No long reads provided, skipping long reads assembly" > LR_assembly.fasta
        echo "GENERA info: No long reads provided, skipping long reads assembly" >> GENERA-Assembler.log
        """
    }
}


//polishing
process polishing {
	//informations

	//input output
    input:
    file 'SPADES_scaffolds.fasta' from spades_ch1
    file "LR_assembly.fasta" from lr_ch1
    file "R1-fastp.fastq" from fastpR1_ch2
    file "R2-fastp.fastq" from fastpR2_ch2
    file "GENERA-Assembler.log" from log_ch3
    val cpu from params.cpu
    val ram from params.RAM
    
    output: 
    file "assembly_final.fasta" into pilon_ch2
    file "assembly_final.fasta" into pilon_ch3
    file "assembly_final.fasta" into pilon_ch4
    file "assembly_final.fasta" into pilon_ch5
    file "GENERA-Assembler.log" into log_ch4


    //script
    script:
    if (params.ontreads != 'no'){
        println "GENERA info: polishing ONT assembly"
        """
        bwa index LR_assembly.fasta
        bwa mem -t $cpu LR_assembly.fasta R1-fastp.fastq R2-fastp.fastq > short_read_mapping.sam
        samtools sort short_read_mapping.sam -@ $cpu -o short_read_mapping.bam
        samtools index short_read_mapping.bam
        java -Xmx500G -jar /opt/pilon-1.24.jar --genome LR_assembly.fasta --bam short_read_mapping.bam --outdir CORR --output assembly_pilon
        mv CORR/assembly_pilon.fasta assembly_final.fasta
        echo "GENERA info: polishing ONT assembly" >> GENERA-Assembler.log
        """
    }
    else if (params.pbreads != 'no'){
        println "GENERA info: polishing Pacbio assembly"
        """
        bwa index LR_assembly.fasta
        bwa mem -t $cpu LR_assembly.fasta R1-fastp.fastq R2-fastp.fastq > short_read_mapping.sam
        samtools sort short_read_mapping.sam -@ $cpu -o short_read_mapping.bam
        samtools index short_read_mapping.bam
        java -Xmx500G -jar /opt/pilon-1.24.jar --genome LR_assembly.fasta --bam short_read_mapping.bam --outdir CORR --output assembly_pilon
        mv CORR/assembly_pilon.fasta assembly_final.fasta
        echo "GENERA info: polishing PB assembly" >> GENERA-Assembler.log
        """
    }
    else {
        println "GENERA info: No Long reads assembly detected, skipping polishing"
        """
        mv SPADES_scaffolds.fasta assembly_final.fasta
        echo "GENERA info: No Long reads assembly detected, skipping polishing" >> GENERA-Assembler.log
        """       
    }
}

//Mapping reads - short reads
process shortreadsmapping {
	//informations

	//input output
    input:
    file "assembly_final.fasta" from pilon_ch2
    file "R1-fastp.fastq" from fastpR1_ch3
    file "R2-fastp.fastq" from fastpR2_ch3
    file "GENERA-Assembler.log" from log_ch4
    val cpu from params.cpu
    
    output:
    file "depth.txt" into depthSR_ch1
    file "depth.txt" into depthSR_ch2
    file "R1-sort.bam" into r1SRbam_ch1
    file "R2-sort.bam" into r2SRbam_ch1
    file "R1-sort.bam.bai" into r1SRbamindex_ch1
    file "R2-sort.bam.bai" into r2SRbamindex_ch1
    file "GENERA-Assembler.log" into log_ch5

    //script
    script:
    """
    bwa index assembly_final.fasta
    #R1
    bwa mem -t $cpu assembly_final.fasta R1-fastp.fastq > R1.sam
    samtools sort R1.sam -@ $cpu -o R1-sort.bam
    samtools index R1-sort.bam
    #R2
    bwa mem -t $cpu assembly_final.fasta R2-fastp.fastq > R2.sam
    samtools sort R2.sam -@ $cpu -o R2-sort.bam
    samtools index R2-sort.bam
    jgi_summarize_bam_contig_depths --outputDepth depth.txt *sort.bam
    echo "GENERA info: Short reads mapping with bwa mem and samtools" >> GENERA-Assembler.log
    """
}

//Mapping reads - long reads
process longreadsmapping {
	//informations

	//input output
    input:
    file "assembly_final.fasta" from pilon_ch3
    file 'ont.fastq' from ontreads_ch2
    file 'pb.fastq' from pbreads_ch2
    file "GENERA-Assembler.log" from log_ch5
    val cpu from params.cpu
    
    output:
    file "sam.cov" into depthLR_ch
    file "GENERA-Assembler.log" into log_ch6

    //script
    script:
    if (params.ontreads != 'no'){
        println "GENERA info: ONT reads mapping to ref"
        """
        bwa index assembly_final.fasta
        bwa mem -x ont2d -t $cpu assembly_final.fasta ont.fastq > ont.sam
        samtools sort ont.sam -@ $cpu -o ont-sort.bam
        samtools coverage ont-sort.bam > sam.cov
        echo "GENERA info: ONT reads mapping with bwa mem and samtools" >> GENERA-Assembler.log
        """  
    }
    else if (params.pbreads != 'no'){
        println "GENERA info: Pacbio reads mapping to ref"
        """
        bwa index assembly_final.fasta
        bwa mem -x pacbio -t $cpu assembly_final.fasta pb.fastq > pb.sam
        samtools sort pb.sam -@ $cpu -o pb-sort.bam
        samtools coverage pb-sort.bam > sam.cov
        echo "GENERA info: PB reads mapping with bwa mem and samtools" >> GENERA-Assembler.log
        """  
    }
    else {
        println "GENERA info: No long reads, skipping long reads mapping to ref"
        """
        echo "GENERA info: No long reads, skipping long reads mapping to ref" > sam.cov
        echo "GENERA info: No long reads, skipping long reads mapping to ref" >> GENERA-Assembler.log
        """  
    }
}

//Binning
process binning {
	//informations

	//input output
    input:
    file "assembly_final.fasta" from pilon_ch4
    file "depth.txt" from depthSR_ch2
    file "R1-sort.bam" from r1SRbam_ch1
    file "R2-sort.bam" from r2SRbam_ch1
    file "R1-sort.bam.bai" from r1SRbamindex_ch1
    file "R2-sort.bam.bai" from r2SRbamindex_ch1
    file "GENERA-Assembler.log" from log_ch6
    val cpu from params.cpu
    
    output:
    file 'METABAT_bin-*.fa' into binsM_ch1
    file 'CONCOCT_bin-*.fa' into binsC_ch1
    file 'METABAT_bin-*.fa' into binsM_ch2
    file 'CONCOCT_bin-*.fa' into binsC_ch2
    file "GENERA-Assembler.log" into log_ch7


    //script
    script:
    if (params.binner == 'metabat'){
        println "GENERA info: Metabat binning activated"
        """
        runMetaBat.sh -t $cpu assembly_final.fasta R1-sort.bam R2-sort.bam
        cd assembly_final.fasta.metabat-bins*/
        find *.fa | cut -f2 -d"." > fa.list
        for f in `cat fa.list`; do cp bin.\$f.fa METABAT_bin-\$f.fa; done
        mv METABAT_bin*.fa ../
        echo "GENERA info: no binning activated" > CONCOCT_bin-FALSE.fa
        echo "GENERA info: Metabat binning activated" >> GENERA-Assembler.log
        """ 
    }
    else if (params.binner == 'concoct'){
        println "GENERA info: CONCOCT binning activated"
        """
        cut_up_fasta.py assembly_final.fasta -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
        concoct_coverage_table.py contigs_10K.bed *sort.bam > coverage_table.tsv
        concoct -t $cpu --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/
        merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
        mkdir concoct_output/fasta_bins
        extract_fasta_bins.py assembly_final.fasta concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
        cd concoct_output/fasta_bins
        find *.fa | cut -f1 -d"." > fa.list
        for f in `cat fa.list`; do cp \$f.fa CONCOCT_bin-\$f.fa; done
        mv CONCOCT_bin* ../../
        echo "GENERA info: no binning activated" > METABAT_bin-FALSE.fa
        echo "GENERA info: CONCOCT binning activated" >> GENERA-Assembler.log
        """ 
    }
    else if (params.binner == 'all') {
        println "GENERA info: Metabat and CONCOCT binning activated"
        """
        #metabat
        runMetaBat.sh -t $cpu assembly_final.fasta R1-sort.bam R2-sort.bam
        cd assembly_final.fasta.metabat-bins*/
        find *.fa | cut -f2 -d"." > fa.list
        for f in `cat fa.list`; do cp bin.\$f.fa METABAT_bin-\$f.fa; done
        mv METABAT_bin*.fa ../
        cd ../
        #concoct
        cut_up_fasta.py assembly_final.fasta -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
        concoct_coverage_table.py contigs_10K.bed *sort.bam > coverage_table.tsv
        concoct -t $cpu --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/
        merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
        mkdir concoct_output/fasta_bins
        extract_fasta_bins.py assembly_final.fasta concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
        cd concoct_output/fasta_bins
        find *.fa | cut -f1 -d"." > fa.list
        for f in `cat fa.list`; do cp \$f.fa CONCOCT_bin-\$f.fa; done
        mv CONCOCT_bin* ../../ 
        echo "GENERA info: Metabat and CONCOCT binning activated" >> GENERA-Assembler.log    
        """ 
    }
    else{
       println "GENERA info: no binning activated" 
        """
        echo "GENERA info: no binning activated" > METABAT_bin-FALSE.fa
        echo "GENERA info: no binning activated" > CONCOCT_bin-FALSE.fa
        echo "GENERA info: no binning activated" >> GENERA-Assembler.log
        """ 
    }
}

//RagTag mapping
process ragtag {
	//informations

	//input output
    input:
    file "assembly_final.fasta" from pilon_ch5
    file 'reference.fasta' from refgenome_ch
    file "GENERA-Assembler.log" from log_ch7
    val cpu from params.cpu
    
    output:
    file "assembly_final.fasta" into ragoo_ch1
    file "assembly_final.fasta" into ragoo_ch2
    file "assembly_final.fasta" into ragoo_ch3
    file "GENERA-Assembler.log" into log_ch8

    //script
    script:
    if (params.ragtag == 'yes'){
        if (params.metagenome == 'yes') {
            println "GENERA info: ragtag options is not authorized with the metagenomic mode"
            """
            """
        }
        if (params.refgenome == 'no') {
            println "GENERA info: a reference genome must be provided with ragtag option"
            """
            """
        }
        else {
            """
            #ref first and genome in second
            ragtag.py scaffold reference.fasta assembly_final.fasta  -t $cpu
            mv ragtag_output/ragtag.scaffold.fasta .
            mv ragtag.scaffold.fasta assembly_final.fasta
            echo "GENERA info: ragtag activated: scaffold mode" >> GENERA-Assembler.log
            """
        }
    }
    else {
        println "GENERA info: ragtag not activated"
        """
        mv assembly_final.fasta temp.fasta
        mv temp.fasta assembly_final.fasta
        echo "GENERA info: ragtag not activated" >> GENERA-Assembler.log
        """   
    }
}

//Quast
process quast {
	//informations
  
	//input output
    input:
    file "assembly_final.fasta" from ragoo_ch3
    file 'METABAT_bin-*.fa' from binsM_ch1
    file 'CONCOCT_bin-*.fa' from binsC_ch1
    file "GENERA-Assembler.log" from log_ch8
    val cpu from params.cpu

    output:
    file "all-quast.results" into quast_ch
    file "GENERA-Assembler.log" into log_ch9


    //script
    script:
    if (params.binner == 'none'){
        println "GENERA info: no binning activated, quast will run on genome"
        """
        quast.py assembly_final.fasta -t $cpu
        cat quast_results/*/report.txt > all-quast.results
        echo "GENERA info: no binning activated, quast will run on genome" >> GENERA-Assembler.log
        """
    }
    else if (params.binner == 'metabat'){
        println "GENERA info: metabat binning activated, quast will run on bins"
        """
        for f in METABAT_bin-*.fa; do quast.py \$f -t $cpu; done
        cat quast_results/*/report.txt > all-quast.results
        echo "GENERA info: metabat binning activated, quast will run on bins" >> GENERA-Assembler.log
        """
    }
    else if (params.binner == 'concoct'){
        println "GENERA info: concoct binning activated, quast will run on bins"
        """
        for f in CONCOCT_bin-*.fa; do quast.py \$f -t $cpu; done
        cat quast_results/*/report.txt > all-quast.results
        echo "GENERA info: concoct binning activated, quast will run on bins" >> GENERA-Assembler.log
        """
    }
    else if (params.binner == 'all'){
        println "GENERA info: all binning activated, quast will run on bins"
        """
        for f in METABAT_bin-*.fa; do quast.py \$f -t $cpu; done
        for f in CONCOCT_bin-*.fa; do quast.py \$f -t $cpu; done
        cat quast_results/*/report.txt > all-quast.results
        echo "GENERA info: all binning activated, quast will run on bins" >> GENERA-Assembler.log
        """
    }

}

//output the results
process publicationResults {
	//informations
    publishDir "$params.outdir", mode: 'copy', overwrite: false

	//input output
    input:
    file 'R1-fastp_fastqc.html' from fastqcR1_ch
    file 'R2-fastp_fastqc.html' from fastqcR2_ch   
    file "assembly_final.fasta" from ragoo_ch1
    file "sam.cov" from depthLR_ch
    file "depth.txt" from depthSR_ch1
    file 'METABAT_bin-*.fa' from binsM_ch2
    file 'CONCOCT_bin-*.fa' from binsC_ch2
    file "all-quast.results" from quast_ch
    file "GENERA-Assembler.log" from log_ch9
    val binner from params.binner
    val companion from params.companion
    val version from params.version


    output:
    file "FASTQC-R1.html" into fastqcR1FINAL_ch
    file "FASTQC-R2.html" into fastqcR2FINAL_ch
    file "Genome.fasta" into genomeFINAL_ch
    file "LongReads-coverage.txt" into depthLRFINAL_ch
    file "ShortReads-coverage.txt" into depthSRFINAL_ch
    file 'METABAT_bin-*.fasta' into binsMFINAL_ch
    file 'CONCOCT_bin-*.fasta' into binsCFINAL_ch
    file "GENERA-Assembler.log" into logFINAL_ch
    file "quast.results" into quastFINAL_ch

    //script
    script:
    if (params.binner != 'none'){
        """
        mv R1-fastp_fastqc.html FASTQC-R1.html
        mv R2-fastp_fastqc.html FASTQC-R2.html
        mv assembly_final.fasta Genome.fasta
        mv sam.cov LongReads-coverage.txt
        mv depth.txt ShortReads-coverage.txt
        mv all-quast.results quast.results
        #Compute percentage of binning
        $companion Genome.fasta --mode=$binner
        mv binned.info GENERA-Assembler.log
        echo "GENERA info: binning activated, binned percentage computed" >> GENERA-Assembler.log
        echo VERSION: >> GENERA-Assembler.log
        echo $version >> GENERA-Assembler.log
        #rename of overwrite
        find CONCOCT_bin-*.fa | cut -f1 -d"." > concoct.list
        for f in `cat concoct.list`; do mv \$f.fa \$f.fasta; done
        find METABAT_bin-*.fa | cut -f1 -d"." > metabat.list
        for f in `cat metabat.list`; do mv \$f.fa \$f.fasta; done
        """
    }
    else{
        """
        mv R1-fastp_fastqc.html FASTQC-R1.html
        mv R2-fastp_fastqc.html FASTQC-R2.html
        mv assembly_final.fasta Genome.fasta
        mv sam.cov LongReads-coverage.txt
        mv depth.txt ShortReads-coverage.txt
        mv all-quast.results quast.results
        echo "GENERA info: no binning activated, no binned percentage computed" >> GENERA-Assembler.log
        echo VERSION: >> GENERA-Assembler.log
        echo $version >> GENERA-Assembler.log
        find CONCOCT_bin-*.fa | cut -f1 -d"." > concoct.list
        for f in `cat concoct.list`; do mv \$f.fa \$f.fasta; done
        find METABAT_bin-*.fa | cut -f1 -d"." > metabat.list
        for f in `cat metabat.list`; do mv \$f.fa \$f.fasta; done
        """
    }
}
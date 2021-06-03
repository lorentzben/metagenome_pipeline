#!/usr/bin/env nextflow

def helpMessage(){
	log.info"""
	Usage:

    This script should be called as follows:
    srun --pty  -p inter_p  --mem=40G --nodes=2 --ntasks-per-node=8 --time=12:00:00 --job-name=qlogin /bin/bash -l
    TODO update this command 
    nextflow run main.nf 

    Main arguments:
        --input [path/to/folder]      Folder containing reads from this analysis
        --outdir [file]               The output directory where the results will be saved 

    """.stripIndent()
}


if(params.input){
    Channel
        .fromPath(params.input)
        .ifEmpty {exit 1, log.info "Cannot find path file ${input}"}
        .set{ ch_seqs }
}

//show help message 
if (params.help){
    helpMessage()
    exit 0 
}

process QualityTrim{

    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    container "docker://lorentzb/trim-galore"

    input:
    path seq_dir from ch_seqs

    output:
    file "done.txt" into ch_done
    path "clean_reads" into ch_cleaned
    

    shell:
    '''
    #!/usr/bin/env bash

    eval "$(conda shell.bash hook)"
    conda activate metagenome  

    trim_galore -q 20 --fastqc --cores 4 -o clean_reads !{seq_dir}/*.fastq

    echo "done" > done.txt

    '''
}

process multiQC{
    publishDir "${params.outdir}/mulit_qc", mode: 'copy'

    container "docker://lorentzb/multi_qc"

    input:
    path "clean_reads" from ch_cleaned

    output:
    file "multiqc_report.html" into ch_multi_qc_report

    script:
    """
    #!/usr/bin/env bash

    multiqc clean_reads/    
    
    """
    
}

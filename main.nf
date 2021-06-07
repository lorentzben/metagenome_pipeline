#!/usr/bin/env nextflow

def helpMessage(){
	log.info"""
	Usage:

    This script should be called as follows:
    srun --pty  -p inter_p  --mem=40G --nodes=2 --ntasks-per-node=8 --time=12:00:00 --job-name=qlogin /bin/bash -l
    TODO update this command 
    nextflow run main.nf --input seqs --mapping mapping.tsv --outdir paired 
    This program expects paired end reads, a single end read pipeline can be developed if required. 

    Main arguments:
        --input [path/to/folder]      Folder containing reads from this analysis
        --mapping [file]              Mapping file with, sample-id \t forward-read \t reverse-read
        --outdir [file]               The output directory where the results will be saved 

    """.stripIndent()
}


if(params.input){
    Channel
        .fromPath(params.input)
        .ifEmpty {exit 1, log.info "Cannot find path file ${input}"}
        .set{ ch_seqs }
}

if(params.mapping){
    Channel
        .fromPath(params.mapping)
        .ifEmpty {exit 1, log.info "Cannot find path file ${mapping}"}
        .into{ ch_mapping_file ; ch_mapping_file_sam }
}

Channel
    .fromPath("${baseDir}/Gallus_gallus")
    .ifEmpty{ exit 1, log.info "Cannot find ref database, please collect"}
    .set{ ch_chicken_ref }

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
    path "clean_reads" into ch_cleaned_contam
    

    shell:
    '''
    #!/usr/bin/env bash

    eval "$(conda shell.bash hook)"
    conda activate metagenome  

    trim_galore -q 20 --fastqc --paired --cores 4 -o clean_reads !{seq_dir}/*.fastq

    echo "done" > done.txt

    '''
}

process multiQC{
    publishDir "${params.outdir}/multi_qc", mode: 'copy'

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

process CheckForContamination{
    publishDir "${params.outdir}/decontaminate", mode: 'copy'

    container "docker://lorentzb/bowtie"

    input:
    file mapping from ch_mapping_file
    path "clean_reads" from ch_cleaned_contam
    path "Gallus_gallus" from ch_chicken_ref

    output:
    path "decontam" into ch_decontam

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    subprocess.run(['mkdir decontam'],shell=True)
    import pandas as pd
    samples = pd.read_csv("${mapping}",sep='\t')

    for index, row in samples.iterrows():
        forward = row['forward-read']
        reverse = row['reverse-read']
        stub = row['sequence-id']
        command = "bowtie2 -p 4 -x Gallus_gallus/genome -1 clean_reads/" +forward[:-6]+"_val_1.fq -2 clean_reads/"+reverse[:-6]+ "_val_2.fq -S decontam/"+stub+".sam"
        result = subprocess.run([command], shell=True)
    
  
    """
}

process ConvertSamtoBamandModify{
    publishDir "${params.outdir}/decontaminate", mode: 'copy'

    container "docker://"

    input:
    file mapping from ch_mapping_file_sam
    path "decontam" from ch_decontam
    
    output:
    path "no_host" into ch_no_host_seqs

    script:
    """
    #!/usr/bin/env python3
    import glob
    import subprocess

    sample_list = glob.glob('./decontam/*.sam')
    subprocess.run(['mkdir no_host'],shell=True)

    for sample in sample_list:
        sample = sample.split('/')[2]
        #convert sam to bam 
        sam_conv = 'samtools view -bS decontam/'+sample+' > '+sample[:-4]+'.bam'
        subprocess.run([sam_conv], shell=True)

        #filter the unmapped reads
        filter_command = 'samtools view -b -f 12 -F 256 '+sample[:-4]+'.bam > '+sample[:-4]+'_bothUnmapped.bam'
        subprocess.run([filter_command],shell=True)

        #split reads into fastq files
        sort_command = 'samtools sort -n -m 5G -@ 2 '+sample[:-4]+'_bothUnmapped.bam -o '+sample[:-4]+'_sorted.bam'
        subprocess.run([sort_command], shell=True)
        split_command = 'samtools fastq -@ 8 '+sample[:-4]+'_sorted.bam -1 no_host/'+sample[:-4]+'_R1.fastq.gz -2 no_host/'+sample[:-4]+'_R2.fastq.gz -0 /dev/null -s /dev/null -n'
        subprocess.run([split_command],shell=True)


    """
}
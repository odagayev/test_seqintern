nextflow.enable.dsl=2

params.reads = 'data/reads/*.fastq'
params.index = 'data/index'
params.gtf = 'data/annotation.gtf'

process fastqc {
    input:
    path reads from params.reads

    output:
    path 'fastqc_reports' into fastqc_results

    script:
    """
    mkdir -p fastqc_reports
    fastqc -o fastqc_reports $reads
    """
}

process hisat2_align {
    input:
    path reads from params.reads

    output:
    path 'aligned.bam' into aligned_bam

    script:
    """
    hisat2 -x ${params.index} -U $reads -S aligned.sam
    samtools view -bS aligned.sam > aligned.bam
    """
}

process featureCounts {
    input:
    path bam from aligned_bam

    output:
    path 'counts.txt' into counts

    script:
    """
    featureCounts -a ${params.gtf} -o counts.txt $bam
    """
}

workflow {
    fastqc()
    hisat2_align()
    featureCounts()
}

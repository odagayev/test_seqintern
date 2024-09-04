nextflow.enable.dsl=2

params.reads = 'data/reads/*.fastq.gz'
params.index = 'data/index'
params.gtf = 'data/annotation.gtf'

process fastqc {
    input:
    path reads from params.reads

    output:
    path 'fastqc_results'

    script:
    """
    fastqc $reads -o fastqc_results
    """
}

process trimmomatic {
    input:
    path reads from params.reads

    output:
    path 'trimmed_reads'

    script:
    """
    trimmomatic PE -phred33 $reads trimmed_reads
    """
}

process hisat2 {
    input:
    path reads from trimmomatic.out

    output:
    path 'aligned_reads'

    script:
    """
    hisat2 -x ${params.index} -U $reads -S aligned_reads
    """
}

process featureCounts {
    input:
    path bam from hisat2.out

    output:
    path 'counts.txt'

    script:
    """
    featureCounts -a ${params.gtf} -o counts.txt $bam
    """
}
}

workflow {
    reads = Channel.fromPath(params.reads)
    fastqc(reads)
    trimmed_reads = trimmomatic(reads)
    aligned_reads = hisat2(trimmed_reads)
    featureCounts(aligned_reads)
}
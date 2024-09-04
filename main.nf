nextflow.enable.dsl=2

params.reads = 'data/reads/*.fastq'
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

process hisat2_align {
    input:
    path reads from params.reads
    path index from params.index

    output:
    path 'aligned.bam'

    script:
    """
    hisat2 -x $index -U $reads -S aligned.sam
    samtools view -bS aligned.sam > aligned.bam
    """
}

process featurecounts {
    input:
    path bam from hisat2_align.out
    path gtf from params.gtf

    output:
    path 'counts.txt'

    script:
    """
    featureCounts -a $gtf -o counts.txt $bam
    """
}

workflow {
    fastqc()
    hisat2_align()
    featurecounts()
}

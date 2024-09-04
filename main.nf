nextflow.enable.dsl=2

params.reads = 'data/reads/*.fastq.gz'
params.genome = 'data/genome/genome.fa'
params.gtf = 'data/genome/annotations.gtf'

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

process trim {
    input:
    path reads from params.reads

    output:
    path 'trimmed_reads'

    script:
    """
    trimmomatic SE -phred33 $reads trimmed_reads ILLUMINACLIP:TruSeq3-SE.fa:2:30:10
    """
}

process align {
    input:
    path reads from trim.out

    output:
    path 'aligned_reads'

    script:
    """
    hisat2 -x ${params.genome} -U $reads -S aligned_reads
    """
}

process quantify {
    input:
    path bam from align.out

    output:
    path 'counts.txt'

    script:
    """
    featureCounts -a ${params.gtf} -o counts.txt $bam
    """
}

workflow {
    reads_ch = Channel.fromPath(params.reads)
    genome_ch = Channel.value(params.genome)
    gtf_ch = Channel.value(params.gtf)

    reads_ch | fastqc
    reads_ch | trim | align | quantify
}

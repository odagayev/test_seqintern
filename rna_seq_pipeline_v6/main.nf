nextflow.enable.dsl=2

params.case_samples = 'path/to/case_samples/*'
params.control_samples = 'path/to/control_samples/*'
params.reference_genome = 'path/to/reference/genome'
params.gtf_file = 'path/to/annotation.gtf'

process fastqc {
    input:
    path sample

    output:
    path 'fastqc_reports'

    script:
    """
    fastqc $sample -o fastqc_reports
    """
}

process alignReads {
    input:
    path sample

    output:
    path 'aligned_reads'

    script:
    """
    hisat2 -x ${params.reference_genome} -U $sample -S aligned_reads/${sample.baseName}.sam
    """
}

process quantifyExpression {
    input:
    path sam_file

    output:
    path 'expression_counts'

    script:
    """
    featureCounts -a ${params.gtf_file} -o expression_counts/${sam_file.baseName}.counts $sam_file
    """
}

process differentialExpression {
    input:
    path case_counts, control_counts

    output:
    path 'differential_expression_results'

    script:
    """
    Rscript -e "
    library(DESeq2)
    case_counts <- read.table('$case_counts', header=TRUE, row.names=1)
    control_counts <- read.table('$control_counts', header=TRUE, row.names=1)
    coldata <- data.frame(condition=factor(c('case', 'control')))
    dds <- DESeqDataSetFromMatrix(countData=cbind(case_counts, control_counts), colData=coldata, design=~condition)
    dds <- DESeq(dds)
    res <- results(dds)
    write.csv(as.data.frame(res), file='differential_expression_results/deseq2_results.csv')
    "
    """
}

workflow {
    Channel
        .fromPath(params.case_samples)
        .set { case_samples }

    Channel
        .fromPath(params.control_samples)
        .set { control_samples }

    case_samples
        | fastqc
        | alignReads
        | quantifyExpression
        | set { case_counts }

    control_samples
        | fastqc
        | alignReads
        | quantifyExpression
        | set { control_counts }

    differentialExpression(case_counts, control_counts)
}
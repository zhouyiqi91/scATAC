process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.6--py39hf95cd2a_1' :
        'biocontainers/cutadapt:4.6--py39hf95cd2a_1' }"

    input:
    tuple val(meta), path(reads)
    val(suffix)

    output:
    tuple val(meta), path("*.trim.fastq${suffix}"), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq${suffix}" : "-o ${prefix}_1.trim.fastq${suffix} -p ${prefix}_2.trim.fastq${suffix}"
    """
    cutadapt \\
        -Z \\
        --cores $task.cpus \\
        $args \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def suffix = gzip ? '.gz': ''
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq${suffix}" : "-o ${prefix}_1.trim.fastq${suffix} -p ${prefix}_2.trim.fastq${suffix}"
    """
    touch ${prefix}.cutadapt.log
    touch ${trimmed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}

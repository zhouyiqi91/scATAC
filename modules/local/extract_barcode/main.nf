process EXTRACT_BARCODE {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta) , path(reads)

    output:
    tuple val(meta), path("${meta.id}*.{fastq, fastq.gz}"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    extract_barcode.py --fq1 ${reads[0]} --fq2 ${reads[1]} --fq3 ${reads[2]} --sample ${meta.id}
    """

    stub:
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension = (args2 ==~ extension_pattern) ? (args2 =~ extension_pattern)[0][2].toLowerCase() : "bam"
    def create_unmapped = ""
    if (meta.single_end) {
        create_unmapped = save_unaligned ? "touch ${prefix}.unmapped.fastq.gz" : ""
    } else {
        create_unmapped = save_unaligned ? "touch ${prefix}.unmapped_1.fastq.gz && touch ${prefix}.unmapped_2.fastq.gz" : ""
    }

    """
    touch ${prefix}.${extension}
    touch ${prefix}.bowtie2.log
    ${create_unmapped}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

}

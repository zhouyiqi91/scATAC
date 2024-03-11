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
    extract_barcode.py --fq1 ${reads[0]} --fq2 ${reads[1]} --fq3 ${reads[2]} --sample ${meta.id} --len_bc ${params.len_bc}
    """

}

process SNAPATAC2 {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta) , path(aligned)

    output:
    tuple val(meta), path("${meta.id}.fragments.gz"), emit: fragments

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    snap.pp.make_fragment_file("${aligned}", "${meta.id}.fragments.gz", barcode_regex="(.*):.*:.*")
    genome = snap.genome.Genome(fasta="${params.fasta}", annotation="${params.gtf}")
    data = snap.pp.import_data(
        "${meta.id}.fragments.gz",
        chrom_sizes=genome,
        sorted_by_barcode=False,
    )

    """

}

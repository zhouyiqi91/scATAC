process SNAPATAC2 {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta) , path(aligned)

    output:
    tuple val(meta), path("${meta.id}.fragments.gz"), emit: fragment
    path('*.png'), emit: image
    path('*.h5ad'), emit: h5ad
    path('*.log'), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    snap.pp.make_fragment_file("${aligned}", "${meta.id}.fragments.gz", barcode_regex="(.*):.*:.*")
    genome = snap.genome.Genome(fasta="${params.fasta}", annotation="${params.gtf}")
    data = snap.pp.import_data(
        "${meta.id}.fragments.gz",
        chrom_sizes=genome,
        file="${meta.id}.snap.h5ad",
        sorted_by_barcode=False,
    )
    snap.pl.frag_size_distr(data, interactive=False, outfile="${meta.id}.fragSize.png")
    snap.metrics.tsse(data, genome)
    snap.pl.tsse(data, interactive=False, outfile="${meta.id}.tsse.png")
    median_tsse = np.median(data.obs['tsse'])
    snap.pp.filter_cells(data, min_counts=1, min_tsse=10, max_counts=100000)
    """

}

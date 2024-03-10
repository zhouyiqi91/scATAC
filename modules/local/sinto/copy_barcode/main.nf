    //TODO
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta) , path(reads)

    output:
    tuple val(meta), path("*.trim.fastq${suffix}"), emit: reads
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def unaligned = ""
    def reads_args = ""
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-U ${reads}"
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
    }

    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension = (args2 ==~ extension_pattern) ? (args2 =~ extension_pattern)[0][2].toLowerCase() : "bam"

    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 \\
        -x \$INDEX \\
        $reads_args \\
        --threads $task.cpus \\
        $unaligned \\
        $args \\
        2> >(tee ${prefix}.bowtie2.log >&2) \\
        | samtools $samtools_command $args2 --threads $task.cpus -o ${prefix}.${extension} -

    if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
        mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
    fi

    if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
        mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """


process EXTRACTSEQNAMES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(seqtsv), path(fasta)

    output:
    tuple val(meta), path("*.seqnames"), emit: seqnames
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zgrep -v '^seqid' ${seqtsv} | \\
        cut -f 1 | \\
        sort | \\
        join -t\$'\\t' - <(grep '>' ${fasta} | sed 's/>//' | sed 's/\\([^ ]\\+\\) .*/\\1\\t&/' | sort) | \\
        cut -f 2 > ${prefix}.correct.seqnames

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        join: \$(join --version | grep '^join' | sed 's/join (GNU coreutils) //')
    END_VERSIONS
    """
}

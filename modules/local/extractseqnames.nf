
process EXTRACTSEQNAMES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.2.3.dev0--0':
        'biocontainers/sed:4.2.3.dev0--0' }"

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
    grep -f <(gunzip -c $seqtsv | cut -f 1) $fasta | sed 's/^>//' > ${prefix}.correct.seqnames

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU sed: \$(echo \$(sed --version | grep '^sed') | sed 's/sed (GNU sed)//')
        GNU grep: \$(grep --version | sed 's/grep (GNU grep) //')
    END_VERSIONS
    """
}

process GTDBFIXNAMES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0':
        'biocontainers/sed:4.7.0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gtdbclean.fna"), emit: fasta
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cat    = fasta.name.endsWith('gz') ? "gunzip -c ${fasta}" : "cat ${fasta}"
    """
    $cat | sed '/>/s/~.*//' > ${prefix}.gtdbclean.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version | grep '^sed') | sed 's/sed (GNU sed)//')
    END_VERSIONS
    """
}

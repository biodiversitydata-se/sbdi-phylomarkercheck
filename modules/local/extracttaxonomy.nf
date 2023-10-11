process EXTRACTTAXONOMY {
    tag "$meta.id"
    label 'process_single'

    conda "conda::sed=4.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0':
        'biocontainers/sed:4.7.0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: taxonomy
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep '>' $fasta | sed 's/>\\([^ ]\\+\\) \\([^ ]\\+\\) .*/\\1\\t\\2/' > ${prefix}.taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version | grep '^sed') | sed 's/sed (GNU sed)//')
    END_VERSIONS
    """
}

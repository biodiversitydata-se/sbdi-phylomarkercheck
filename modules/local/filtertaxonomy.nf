
process FILTERTAXONOMY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta), path(taxonomy)

    output:
    tuple val(meta), path(fasta, includeInputs: true), path("*.filttax.tsv"), emit: filtered, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ \$( grep -c '>' $fasta ) -gt 0 ]; then
        grep '>' $fasta | sed 's/>//' | sort | join -t\$'\\t' - <(cat $taxonomy | sort -t\$'\\t' -k 1b,1) > ${prefix}.filttax.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        join: \$(join --version | grep '^join' | sed 's/join (GNU coreutils) //')
    END_VERSIONS
    """
}

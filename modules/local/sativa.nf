
process SATIVA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::sativa=0.9.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quay.io/biocontainers/sativa:0.9.3--py312h031d066_0':
        'biocontainers/sativa:0.9.3--py312h031d066_0'
    }"

    input:
    tuple val(meta), path(alignment), path(taxonomy)

    output:
    tuple val(meta), path("*.mis")             , emit: misplaced
    tuple val(meta), path("*.log")             , emit: log
    tuple val(meta), path("*.refjson")         , emit: refjson
    tuple val(meta), path("*.final_epa.jplace"), emit: final_epa
    tuple val(meta), path("*.l1out_seq.jplace"), emit: l1out_seq
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir tmp
    sativa.py \\
        -s $alignment \\
        -t $taxonomy \\
        -tmpdir tmp \\
        -n ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sativa: 0.9.3
    END_VERSIONS
    """
}

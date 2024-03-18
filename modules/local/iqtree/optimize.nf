process IQTREE_OPTIMIZE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.2.6--h21ec9f0_0':
        'biocontainers/iqtree:2.2.6--h21ec9f0_0' }"

    input:
    tuple val(meta), path(aln), path(tree)

    output:
    tuple val(meta), path("*.brlenopt.treefile")     , emit: phylo
    tuple val(meta), path("*.brlenopt.iqtree")       , emit: iqtree
    tuple val(meta), path("*.brlenopt.log")          , emit: log
    tuple val(meta), path("*.brlenopt.uniqueseq.phy"), emit: uniq, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory.toString().replaceAll(' ', '')
    """
    iqtree \\
        $args \\
        -s $aln \\
        -te $tree \\
        -nt AUTO \\
        -mem $memory \\
        -ntmax $task.cpus \\
        -pre ${prefix}.brlenopt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}

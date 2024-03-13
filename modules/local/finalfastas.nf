process FINALFASTAS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0':
        'biocontainers/sed:4.7.0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.assignTaxonomy.fna.gz"), emit: assign_taxonomy
    tuple val(meta), path("*.addSpecies.fna.gz")    , emit: add_species
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c $fasta | \\
        sed '/^>/s/>[^ ]\\+ \\(d__[^;]\\+\\);\\(.*\\) \\[location.*/>\\1;\\1;\\2/' | \\
        sed '/^>/s/[a-z]__//g' | \\
        gzip -c > ${prefix}.assignTaxonomy.fna.gz

    gunzip -c $fasta | \\
        sed '/^>/s/\\(>[^ ]\\+\\) .*;s__\\(.*\\) \\[location.*/\\1 \\2/' | \\
        gzip -c > ${prefix}.addSpecies.fna.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version | grep '^sed') | sed 's/sed (GNU sed)//')
    END_VERSIONS
    """
}

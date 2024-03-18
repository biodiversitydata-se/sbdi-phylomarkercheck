process BIOPYTHON_REMOVEN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.non.fna"), emit: fasta
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python

    import sys
    import Bio
    from Bio import SeqIO

    output = open("${prefix}.non.fna", "w")
    for record in SeqIO.parse("$fasta", "fasta"):
        if record.seq.count('N') == 0:
            output.write(record.format("fasta"))
    output.close()

    versions = open("versions.yml", "w")
    versions.write(f'"${task.process}":\\n    biopython: {Bio.__version__}\\n')
    versions.close()
    """
}

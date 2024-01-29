
process BIOPYTHON_FILTERGAPPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quay.io/biocontainers/biopython:1.81':
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gapfiltered.fna"), emit: fasta
    path "versions.yml"                     , emit: versions

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

    input    = open("$fasta", "r")
    output   = open("${prefix}.gapfiltered.fna", "w")

    drop_cutoff = 1 - float($args)

    for seqs in SeqIO.parse(input, 'fasta'):
        name = seqs.id
        seq  = seqs.seq

        seqlen = len(seq)
        gapc   = 0

        for i in range(seqlen):
            if seqs[i] == '-':
                gapc += 1

        if ( gapc/float(seqlen) <= drop_cutoff ):
            SeqIO.write(seqs, output, 'fasta')

    input.close()
    output.close()

    versions = open("versions.yml", "w")
    versions.write(f'"${task.process}":\\n    biopython: {Bio.__version__}\\n')
    #versions.write(f"biopython: {Bio.__version__}\\n")
    versions.close()
    """
}

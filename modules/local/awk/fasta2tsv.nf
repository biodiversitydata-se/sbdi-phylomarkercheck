process AWK_FASTA2TSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.1.3--1':
        'biocontainers/gawk:4.1.3--1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: bam
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def read_input = fasta.endsWith('.gz') ? "gunzip -c ${fasta}" : "cat ${fasta}"

    """
    $read_input | \
        awk '/^>/ { printf("\\n%s\\n",\$0); next; } { printf("%s",\$0);} END { printf("\\n");}' | \
        sed '/^\$/d' | paste - -" | \
        gzip -c > ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | grep "GNU Awk" | sed 's/GNU Awk \([0-9.]\+\).*/\1/')
    END_VERSIONS
    """
}

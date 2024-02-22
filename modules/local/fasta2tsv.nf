process FASTA2TSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def read_input = fasta.name.endsWith('.gz') ? "gunzip -c ${fasta}" : "cat ${fasta}"

    """
    $read_input | \
        awk '/^>/ { printf("\\n%s\\n",\$0); next; } { printf("%s",\$0);} END { printf("\\n");}' | \
        sed 's/^>//' | \
        sed '/^\$/d' | \
        paste - - > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | grep "GNU Awk" | sed 's/GNU Awk \\([0-9.]\\+\\).*/\\1/')
    END_VERSIONS
    """
}

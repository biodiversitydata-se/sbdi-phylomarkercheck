process SUBSETTREE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-treeio:1.26.0--r43hdfd78af_1':
        'biocontainers/bioconductor-treeio:1.26.0--r43hdfd78af_1' }"

    input:
    tuple val(meta), path(seqtsv)
    path  phylo

    output:
    tuple val(meta), path("*.subset.newick"), emit: phylo
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(treeio)

    phylo <- read.tree("$phylo")
    names <- read.delim("$seqtsv", sep = '\\t', col.names = c('seqname'), header = FALSE)

    names\$seqname <- sub('~.*', '', names\$seqname)

    ssphylo <- drop.tip(
        phylo,
        phylo\$tip.label[!phylo\$tip.label %in% names\$seqname]
    )

    write.tree(ssphylo, "${prefix}.subset.newick")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    Treeio: ", packageVersion("treeio")) ), "versions.yml")
    """
}

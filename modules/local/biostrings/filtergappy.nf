process BIOSTRINGS_FILTERGAPPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-biostrings:2.70.1--r43ha9d7317_1':
        'biocontainers/bioconductor-biostrings:2.70.1--r43ha9d7317_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gapfiltered.fna"), emit: fasta
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(Biostrings)

    drop_cutoff <- 1 - $args

    seqs <- readRNAStringSet('${fasta}')

    s <- data.frame(
      seqname = names(seqs),
      sequence = seqs
    )
    s\$genome   = sub('~.*', '',  s\$seqname)
    s\$aligned  = gsub('-',   '', s\$sequence)
    s\$naligned = nchar(s\$aligned)

    s <- s[s\$naligned > nchar(s\$sequence) * drop_cutoff,]

    t <- merge(
      s, 
      aggregate(s\$naligned, by = list(s\$genome), FUN = max),
      by.x = c('genome', 'naligned'), by.y = c('Group.1', 'x')
    )[,c('seqname', 'sequence')]

    t\$f <- sprintf(">%s\n%s", t\$seqname, t\$sequence)

    write.table(t\$f, '${prefix}.gapfiltered.fna', quote = FALSE, row.names = FALSE, col.names = FALSE)

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    Biostrings: ", packageVersion("Biostrings")) ), "versions.yml")
    """
}

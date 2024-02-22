process BIOSTRINGS_FILTERGAPPY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-base=4.2.3 conda-forge::r-rmarkdown=2.22 conda-forge::r-tidyverse=2.0.0 conda-forge::r-knitr=1.43 conda-forge::r-dt=0.28 conda-forge::r-dtplyr=1.3.1 conda-forge::r-formattable=0.2.1 conda-forge::r-purrr=1.0.1 conda-forge::r-vegan=2.6_4 conda-forge::r-optparse=1.7.3 conda-forge::r-ggplot2=3.4.2 conda-forge::r-dplyr=1.1.2 conda-forge::r-data.table=1.14.8 conda-forge::r-patchwork=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

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

    def read_input = fasta.endsWith('.gz') ? "gunzip -c ${fasta}" : "cat ${fasta}"
    read_input    += " | awk '/^>/ { printf(\"\\n%s\\n\",\$0); next; } { printf(\"%s\",\$0);} END { printf(\"\\n\");}' | sed '/^\$/d' | paste - -"

    """
    #!/usr/bin/env Rscript

    library(data.table)

    drop_cutoff <- 1 - $args

    s <- fread(cmd = "$read_input", col.names = c('seqname', 'sequence'))

    #s <- data.frame(
    #  seqname = names(seqs),
    #  sequence = seqs
    #)
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

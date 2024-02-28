process UNIQSEQUENCES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-dtplyr=1.3.1 conda-forge::r-data.table=1.14.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.uniq.fna")  , emit: fasta
    tuple val(meta), path("*.uniq.names"), emit: names
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)

    sequences <- read_tsv("${tsv}", col_names = c('name', 'sequence'), col_types = 'cc')

    sequences %>%
        group_by(sequence) %>%
        filter(row_number() == 1) %>%
        ungroup() -> uniq

    uniq %>%
        transmute(s = sprintf(">%s\\n%s", name, sequence)) %>%
        write.table('${prefix}.uniq.fna', quote = FALSE, row.names = FALSE, col.names = FALSE)

    uniq %>%
        select(name) %>%
        write_tsv("${prefix}.uniq.names", col_names = FALSE)

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ",
            paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    tidyr: ", packageVersion('tidyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """
}

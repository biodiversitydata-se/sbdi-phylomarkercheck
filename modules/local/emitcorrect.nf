
process EMITCORRECT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-dtplyr=1.3.1 conda-forge::r-data.table=1.14.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(misplaced)
    path(taxonomy)
    path(metadata)

    output:
    tuple val(meta), path("*.correct.tsv.gz"), emit: correct
    path "versions.yml"                      , emit: versions

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


    misplaced <- read_tsv(c("${misplaced.join('","')}"), skip = 4, col_types = 'ccccdccc') %>%
        rename_with(str_to_lower) %>%
        rename(seqid = `;seqid`) %>%
        distinct(seqid)

    taxonomy <- read_tsv('$taxonomy', col_names = c('seqid', 'taxonomy'), col_types = 'cc') %>%
        mutate(
            t = str_remove_all(taxonomy, '[a-z]__'),
            accession = str_remove(seqid, '~.*')
        ) %>%
        separate(t, c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';')

    metadata <- read_tsv('$metadata', col_types = cols(checkm_completeness = col_double(), checkm_contamination = col_double(), .default = col_character())) %>%
        transmute(accession, quality = checkm_completeness - 5 * checkm_contamination)

    taxonomy %>%
        anti_join(misplaced, by = join_by(seqid)) %>%
        inner_join(metadata, by = join_by(accession)) %>%
        group_by(species) %>%
        slice_max(order_by = quality, n = ${meta.n}) %>%
        write_tsv("${prefix}.correct.tsv.gz")

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

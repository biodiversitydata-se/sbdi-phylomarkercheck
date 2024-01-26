
process EMITCORRECT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-dtplyr=1.3.1 conda-forge::r-data.table=1.14.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(misplacedfwd), path(misplacedrev)
    path(taxonomy)

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

    misplaced <- read_tsv('$misplacedfwd', skip = 4, col_types = 'ccccdccc') %>%
        union(read_tsv('$misplacedrev', skip = 4, col_types = 'ccccdccc')) %>%
        rename_with(str_to_lower) %>%
        rename(seqid = `;seqid`) %>%
        distinct(seqid)

    taxonomy <- read_tsv('$taxonomy', col_names = c('seqid', 'taxonomy'), col_types = 'cc') %>%
        mutate(taxonomy = str_remove_all(taxonomy, '[a-z]__')) %>%
        separate(taxonomy, c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';', remove = FALSE)

    taxonomy %>%
        anti_join(misplaced, by = join_by(seqid)) %>%
        group_by(species) %>%
        slice_sample(n = ${meta.n}) %>%
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

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emitcorrect: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}

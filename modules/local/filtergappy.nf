process FILTERGAPPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1':
        'biocontainers/r-tidyverse:1.2.1' }"

    input:
    tuple val(meta), path(tsv)
    path  metadatatsv

    output:
    tuple val(meta), path("*.gapfiltered.fna"), emit: fasta
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 1.0
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    drop_cutoff <- $args

    s <- read.delim('$tsv', col.names = c('seqname', 'sequence'), sep = '\\t')

    m <- read.delim('$metadatatsv', sep = '\\t') %>%
        transmute(
            genome = accession, gtdb_taxonomy = str_remove_all(gtdb_taxonomy, '[a-z]__'),
            genome_size, checkm_completeness, checkm_contamination, gtdb_representative,
            ncbi_genome_category, ssu_silva_taxonomy
        ) %>%
        separate(gtdb_taxonomy, c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';', extra = 'drop') %>%
        separate(ssu_silva_taxonomy, c('silva_domain', 'silva_phylum', 'silva_class', 'silva_order', 'silva_family', 'silva_genus', 'silva_species'), sep = ';', fill = 'right')

    s %>%
        # Remove sequences with too many gaps
        filter(str_remove_all(sequence, '-') %>% str_length() > str_length(sequence) * drop_cutoff) %>%
        mutate(genome = str_remove(seqname, '~.*')) %>%
        inner_join(m, by = 'genome') %>%
        # Keep the longest sequence, i.e. fewest gaps, from each genome -- 5 for species representative genomes
        group_by(species) %>%
        mutate(
            multiplier = case_when(
                ncbi_genome_category == 'derived from metagenome'  ~ 1,
                ncbi_genome_category == 'derived from single cell' ~ 1,
                gtdb_representative == 't'                         ~ 3,
                TRUE                                               ~ 1
            )
        ) %>%
        mutate(multiplier = ifelse(order != silva_order & silva_order != 'Incertae Sedis' & ! is.na(silva_order), multiplier / 3, multiplier)) %>%
        arrange(str_remove_all(sequence, '-') %>% str_length() %>% desc()) %>%
        filter(row_number() <= 50) %>%
        ungroup() %>%
#        # Select the 50 "best" sequences from each species, preferring GTDB representative genomes and longer sequences
#        mutate(sortnum = ifelse(gtdb_representative == 't', 10, 1) * str_remove_all(sequence, '-') %>% str_length()) %>%
#        group_by(species) %>%
#        arrange(desc(sortnum)) %>%
#        filter(row_number() <= 50) %>%
#        ungroup() %>%
        transmute(s = sprintf(">%s\\n%s", seqname, sequence)) %>%
        write.table('${prefix}.gapfiltered.fna', quote = FALSE, row.names = FALSE, col.names = FALSE)

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    Tidyverse: ", packageVersion("tidyverse")) ), "versions.yml")
    """
}

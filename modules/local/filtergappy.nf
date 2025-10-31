process FILTERGAPPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546'
    }"

    input:
    tuple val(meta), path(tsv)
    path  metadatatsv
    path  ncbi_genome_data

    output:
    tuple val(meta), path("*.gapfiltered.fna"), emit: fasta
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 1.0
    def prefix = task.ext.prefix ?: "${meta.id}"

    def read_ncbi_genome_data = ""
    if ( ncbi_genome_data ) {
        read_ncbi_genome_data =
        """
        ncbi_genome_data = read_tsv(c('${ncbi_genome_data.join("','")}'), skip = 1) %>%
            rename(ncbi_accession = `#assembly_accession`)
        """
    } else {
        read_ncbi_genome_data = "ncbi_genome_data <- tibble(ncbi_accession = character(), excluded_from_refseq = character())"
    }

    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    drop_cutoff <- $args

    s <- read.delim('$tsv', col.names = c('seqname', 'sequence'), sep = '\\t')

    ${read_ncbi_genome_data}

    m <- read.delim('$metadatatsv', sep = '\\t') %>%
        transmute(
            genome = accession, ncbi_accession = str_remove(genome, '^.._'), gtdb_taxonomy = str_remove_all(gtdb_taxonomy, '[a-z]__'),
            genome_size, checkm_completeness, checkm_contamination, gtdb_representative, ncbi_genome_category, ssu_silva_taxonomy
        ) %>%
        anti_join(ncbi_genome_data %>% filter(str_detect(excluded_from_refseq, 'contaminated')), by = join_by(ncbi_accession)) %>%
        separate(gtdb_taxonomy, c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ';', extra = 'drop') %>%
        separate(ssu_silva_taxonomy, c('silva_domain', 'silva_phylum', 'silva_class', 'silva_order', 'silva_family', 'silva_genus', 'silva_species'), sep = ';', fill = 'right', extra = 'drop')

    N_PER_SPECIES = 30
    s %>%
        # Remove sequences with too many gaps
        filter(str_remove_all(sequence, '-') %>% str_length() > str_length(sequence) * drop_cutoff) %>%
        mutate(genome = str_remove(seqname, '~.*')) %>%
        inner_join(m, by = 'genome') %>%
        # We're selecting genes on two criteria: Number per genome (3 for species-representatives, 1 for others) and number per species based on assumed quality (see comments below)
        mutate(
            # ngenes is used to limit the number of sequences to pick from each genome
            ngenes = case_when(
                ncbi_genome_category == 'derived from metagenome'  ~ 1,
                ncbi_genome_category == 'derived from single cell' ~ 1,
                gtdb_representative == 't'                         ~ 3,
                TRUE                                               ~ 1
            ),
            # multiplier is used to weigh each genome within the species
            multiplier = ngenes
        ) %>%
        # Lower the multiplier for mismatches between GTDB order and Silva SSU order, unless the Silva order is NA or 'Incertae Sedis'
        # In cases where the GTDB and Silva taxonomies match, this is fine, in cases where they don't it's affecting the whole species equally and this will have no effect
        mutate(multiplier = ifelse(order != silva_order & silva_order != 'Incertae Sedis' & ! is.na(silva_order), multiplier / 3, multiplier)) %>%
        # Lower the multiplier for more contaminated genomes (even if this is measured with protein coding genes it might say something about the general quality)
        mutate(multiplier = multiplier * checkm_contamination / 100) %>%
        mutate(seq_length = str_remove_all(sequence, '-') %>% str_length()) %>%
        # Select the appropriate number of genes per genome, preferring longer sequences
        group_by(genome) %>%
        arrange(desc(seq_length)) %>%
        filter(row_number() <= ngenes) %>%
        ungroup() %>%
        # Select the top N_PER_SPECIES
        group_by(species) %>%
        arrange(desc(seq_length * multiplier)) %>%
        filter(row_number() <= N_PER_SPECIES) %>%
        ungroup() %>%
#        # Select the 50 "best" sequences from each species, preferring GTDB representative genomes and longer sequences
#        mutate(sortnum = ifelse(gtdb_representative == 't', 10, 1) * str_remove_all(sequence, '-') %>% str_length()) %>%
#        group_by(species) %>%
#        arrange(desc(sortnum)) %>%
#        filter(row_number() <= 50) %>%
#        ungroup() %>%
        transmute(s = sprintf(">%s\\n%s", seqname, sequence)) %>%
        write.table('${prefix}.gapfiltered.fna', quote = FALSE, row.names = FALSE, col.names = FALSE)

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    tidyr: ", packageVersion('tidyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """
}

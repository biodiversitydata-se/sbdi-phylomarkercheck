/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowPhylomarkercheck.initialise(params, log)

// Initialize channels from parameters
Channel.fromPath(params.input)
    .map { [ [ id: "${params.markername}-fwd" ], it ] }
    .set { ch_input }
Channel.fromPath(params.hmm)
    .set { ch_hmm }
Channel.fromPath(params.gtdb_metadata)
    .set { ch_gtdb_metadata }
Channel.of(params.n_per_species.split(','))
    .set { ch_n_per_species }
if ( params.phylogeny ) {
    Channel.fromPath(params.phylogeny)
        .set { ch_phylogeny }
    // Make sure 1 is in ch_n_per_species
    ch_n_per_species
        .concat(Channel.of(1))
        .unique()
}

ch_ncbi_genome_data = Channel.empty()
if ( params.ncbi_genome_data ) {
    ch_ncbi_genome_data = Channel.of(params.ncbi_genome_data.split(','))
        .map { f -> file(f) }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BIOPYTHON_REMOVEN                     } from '../modules/local/biopython/removen'
include { EXTRACTTAXONOMY                       } from '../modules/local/extracttaxonomy'
include { FILTERGAPPY                           } from '../modules/local/filtergappy'
include { FASTA2TSV as ALIGNED2TSV              } from '../modules/local/fasta2tsv'
include { FASTA2TSV as SELECTED2TSV             } from '../modules/local/fasta2tsv'
include { FASTA2TSV as N1SPREP2TSV              } from '../modules/local/fasta2tsv'
include { FILTERTAXONOMY                        } from '../modules/local/filtertaxonomy'
include { SATIVA                                } from '../modules/local/sativa'
include { EMITCORRECT                           } from '../modules/local/emitcorrect.nf'
include { EXTRACTSEQNAMES                       } from '../modules/local/extractseqnames'
include { FINALFASTAS                           } from '../modules/local/finalfastas'
include { UNIQSEQUENCES                         } from '../modules/local/uniqsequences'
include { SUBSETTREE                            } from '../modules/local/subsettree'
include { GTDBFIXNAMES                          } from '../modules/local/gtdbfixnames'
include { IQTREE_OPTIMIZE                       } from '../modules/local/iqtree/optimize'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SEQKIT_SEQ as MAXLEN                  } from '../modules/nf-core/seqkit/seq/main'
include { EMBOSS_REVSEQ                         } from '../modules/nf-core/emboss/revseq/main'
include { HMMER_HMMFETCH                        } from '../modules/nf-core/hmmer/hmmfetch/main'
include { HMMER_HMMALIGN as ALIGNINPUT          } from '../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_HMMALIGN as ALIGNSELECTED       } from '../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_ESLALIMASK as MASKINPUT         } from '../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLALIMASK as MASKSELECTED      } from '../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLREFORMAT as REFORMATINPUT    } from '../modules/nf-core/hmmer/eslreformat/main'
include { HMMER_ESLREFORMAT as REFORMATSELECTED } from '../modules/nf-core/hmmer/eslreformat/main'
include { GUNZIP                                } from '../modules/nf-core/gunzip/main'
include { SEQTK_SUBSEQ as SUBSEQ_SETS           } from '../modules/nf-core/seqtk/subseq/main'
include { SEQTK_SUBSEQ as SUBSEQ_N1SPREPS       } from '../modules/nf-core/seqtk/subseq/main'
include { CAT_CAT                               } from '../modules/nf-core/cat/cat/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PHYLOMARKERCHECK {

    ch_versions = Channel.empty()

    // 0. Remove sequences that are too long or contain Ns
    MAXLEN ( ch_input )
    ch_versions = ch_versions.mix(MAXLEN.out.versions)

    BIOPYTHON_REMOVEN ( MAXLEN.out.fastx )
    ch_versions = ch_versions.mix(BIOPYTHON_REMOVEN.out.versions)

    ch_correct_sequences = BIOPYTHON_REMOVEN.out.fasta

    // 1. Extract the taxonomy from the unaligned fasta and clean up the fasta by removing the taxonomy
    EXTRACTTAXONOMY ( ch_correct_sequences )
    ch_versions = ch_versions.mix(EXTRACTTAXONOMY.out.versions)

    // 2. Reverse the input sequences
    EMBOSS_REVSEQ ( EXTRACTTAXONOMY.out.stripped_fasta.map { [ [ id: "${it[0].id - ~/-fwd/}-rev" ], it[1] ] } )
    ch_versions = ch_versions.mix(EMBOSS_REVSEQ.out.versions)

    // 3. Align both the standard and reverse fasta files with HMMER; filter with rfmask and convert to fasta
    if ( params.hmmkey ) {
        HMMER_HMMFETCH ( ch_hmm.map { [ [ id: params.hmmkey ], it ] }, Channel.of(params.hmmkey), [], [] )
        ch_versions = ch_versions.mix(EXTRACTTAXONOMY.out.versions)
        ch_hmm = HMMER_HMMFETCH.out.hmm.map { it[1] }
    }

    ALIGNINPUT ( EXTRACTTAXONOMY.out.stripped_fasta.mix(EMBOSS_REVSEQ.out.revseq), ch_hmm.first() )
    ch_versions = ch_versions.mix(ALIGNINPUT.out.versions)

    MASKINPUT ( ALIGNINPUT.out.sthlm.map { [ it[0], it[1], [], [], [], [], [], [] ] }, [] )
    ch_versions = ch_versions.mix(MASKINPUT.out.versions)

    REFORMATINPUT(MASKINPUT.out.maskedaln)
    ch_versions = ch_versions.mix(REFORMATINPUT.out.versions)

    ALIGNED2TSV(REFORMATINPUT.out.seqreformated)
    ch_versions = ch_versions.mix(ALIGNED2TSV.out.versions)

    FILTERGAPPY(ALIGNED2TSV.out.tsv, ch_gtdb_metadata.first(), ch_ncbi_genome_data.collect())
    ch_versions = ch_versions.mix(FILTERGAPPY.out.versions)

    // 4. Call sativa with the taxonomy and each filtered alignment
    FILTERGAPPY.out.fasta
        .combine(EXTRACTTAXONOMY.out.taxonomy.map { it[1] })
        .set { ch_gapfiltered_taxonomy }

    FILTERTAXONOMY(ch_gapfiltered_taxonomy)
    ch_versions = ch_versions.mix(FILTERTAXONOMY.out.versions)

    SATIVA(FILTERTAXONOMY.out.filtered)
    ch_versions = ch_versions.mix(SATIVA.out.versions)

    // 5. Emit presumed correct sequences
    ch_n_per_species
        .map { n -> [ id: "${params.markername}-n${n}", n: n ] }
        .combine(SATIVA.out.misplaced.collect { it[1] }.map { [ it ] })
        .set { ch_emitcorrect }

    EMITCORRECT(ch_emitcorrect, FILTERTAXONOMY.out.filtered.map { it[2] }.first(), ch_gtdb_metadata.first())
    ch_versions = ch_versions.mix(EMITCORRECT.out.versions)

    EXTRACTSEQNAMES(
        EMITCORRECT.out.correct
            .combine(ch_correct_sequences.map { it[1] })
    )
    ch_versions = ch_versions.mix(EXTRACTSEQNAMES.out.versions) 

    SUBSEQ_SETS(
        ch_correct_sequences.first().map { it[1] },
        EXTRACTSEQNAMES.out.seqnames
    )
    ch_versions = ch_versions.mix(SUBSEQ_SETS.out.versions)

    FINALFASTAS(SUBSEQ_SETS.out.sequences)
    ch_versions = ch_versions.mix(FINALFASTAS.out.versions)

    // 6. Subset and reoptimize tree
    if ( params.phylogeny ) {
        // Construct a channel with names of sequences for species representatives. Will be used to subset
        // the n1 set for the phylogeny.
        ch_gtdb_metadata
            .splitCsv(header: true, sep: '\t')
            .filter { it.gtdb_representative == 't' }
            .map { [ it.accession ] }
            .join(
                EXTRACTSEQNAMES.out.seqnames
                    .filter { it[0].n == '1' }
                    .map { it[1] }
                    .splitCsv(header: ['seqname'])
                    .map { [ it.seqname - ~/~.*/, it.seqname ] }
            )
            .map { it[1] }
            .collectFile(name: "${params.markername}.n1sprep.txt", newLine: true)
            .map { [ [ id: "${params.markername}-sprep" ], it ] }
            .set { ch_sprep_names }

        // Write a taxonomy file to use in phyloplacement
        ch_sprep_names
            .map { it[1] }
            .splitCsv(header: ['seqname', 'taxonomy', 'sp'], sep: ' ')
            .map { "${it.seqname - ~/~.*/}\t$it.taxonomy $it.sp" }
            .collectFile(name: "${params.markername}-sprep.taxonomy.tsv", newLine: true, storeDir: "${params.outdir}/iqtree")

        // We need a file with non-gappy sequences, i.e. aligned on the correct strand
        CAT_CAT(
            FILTERGAPPY.out.fasta
                .map { it[1] }
                .collect()
                .map{ [ [ id: "${params.markername}.correct" ], it ] }
        )
        ch_versions = ch_versions.mix(CAT_CAT.out.versions)

        SUBSEQ_N1SPREPS(
            CAT_CAT.out.file_out.map { it[1] },
            ch_sprep_names
        )
        ch_versions = ch_versions.mix(SUBSEQ_N1SPREPS.out.versions)

        N1SPREP2TSV(SUBSEQ_N1SPREPS.out.sequences)
        ch_versions = ch_versions.mix(N1SPREP2TSV.out.versions)

        UNIQSEQUENCES(N1SPREP2TSV.out.tsv)
        ch_versions = ch_versions.mix(UNIQSEQUENCES.out.versions)
        
        SUBSETTREE(UNIQSEQUENCES.out.names, ch_phylogeny.first())
        ch_versions = ch_versions.mix(SUBSETTREE.out.versions)

        GTDBFIXNAMES(UNIQSEQUENCES.out.fasta)
        ch_versions = ch_versions.mix(GTDBFIXNAMES.out.versions)

        GTDBFIXNAMES.out.fasta
            .join(SUBSETTREE.out.phylo)
            .set { ch_aln_tree }

        IQTREE_OPTIMIZE(ch_aln_tree)
        ch_versions = ch_versions.mix(IQTREE_OPTIMIZE.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPhylomarkercheck.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPhylomarkercheck.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

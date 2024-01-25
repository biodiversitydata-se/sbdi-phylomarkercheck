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
include { EXTRACTTAXONOMY             } from '../modules/local/extracttaxonomy'
include { SATIVA                      } from '../modules/local/sativa'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { EMBOSS_REVSEQ                 } from '../modules/nf-core/emboss/revseq/main'
include { HMMER_HMMFETCH                } from '../modules/nf-core/hmmer/hmmfetch/main'
include { HMMER_HMMALIGN                } from '../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_ESLALIMASK              } from '../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLREFORMAT             } from '../modules/nf-core/hmmer/eslreformat/main'
include { GUNZIP                        } from '../modules/nf-core/gunzip/main'
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PHYLOMARKERCHECK {

    ch_versions = Channel.empty()

    // 1. Extract the taxonomy from the unaligned fasta and clean up the fasta by removing the taxonomy
    EXTRACTTAXONOMY ( ch_input )
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

    HMMER_HMMALIGN ( EXTRACTTAXONOMY.out.stripped_fasta.mix(EMBOSS_REVSEQ.out.revseq), ch_hmm.first() )
    ch_versions = ch_versions.mix(HMMER_HMMALIGN.out.versions)

    HMMER_ESLALIMASK ( HMMER_HMMALIGN.out.sthlm.map { [ it[0], it[1], [], [], [], [], [], [] ] }, [] )
    ch_versions = ch_versions.mix(HMMER_ESLALIMASK.out.versions)

    HMMER_ESLREFORMAT(HMMER_ESLALIMASK.out.maskedaln)
    ch_versions = ch_versions.mix(HMMER_ESLREFORMAT.out.versions)

    GUNZIP(HMMER_ESLREFORMAT.out.seqreformated)
    ch_versions = ch_versions.mix(GUNZIP.out.versions)

    // 4. Call sativa with the taxonomy and each filtered alignment
    GUNZIP.out.gunzip
        .combine(EXTRACTTAXONOMY.out.taxonomy.map { it[1] })
        .set { ch_sativa }

    SATIVA(ch_sativa)
    ch_versions = ch_versions.mix(SATIVA.out.versions)

    // n.

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

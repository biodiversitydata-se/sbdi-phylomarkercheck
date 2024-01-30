
process SATIVA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::sativa=0.9.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quay.io/biocontainers/sativa:0.9.3--py312h031d066_0':
        'biocontainers/sativa:0.9.3--py312h031d066_0'
    }"

    input:
    tuple val(meta), path(alignment), path(taxonomy)

    output:
    tuple val(meta), path("*.filttax.tsv")     , emit: filtered_taxonomy, optional: true
    tuple val(meta), path("*.mis")             , emit: misplaced
    tuple val(meta), path("*.log")             , emit: log      , optional: true
    tuple val(meta), path("*.refjson")         , emit: refjson  , optional: true
    tuple val(meta), path("*.final_epa.jplace"), emit: final_epa, optional: true
    tuple val(meta), path("*.l1out_seq.jplace"), emit: l1out_seq, optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    if [ \$( grep -c '>' $alignment ) -gt 0 ]; then
        echo "Input not empty, running sativa"
        grep -f <(grep '>' $alignment | sed 's/>//') $taxonomy > ${prefix}.filttax.tsv
        mkdir tmp
        sativa.py \\
            -s $alignment \\
            -t ${prefix}.filttax.tsv \\
            -tmpdir tmp \\
            -n ${prefix} \\
            -T $task.cpus \\
            $args
    else
        cat <<-END_SATIVA > ${prefix}.mis
    ;WARNING: The revised taxon name suggested here is not necessarily the one that has priority in nomenclature. 
    ;Our suggestion should only be taken as indicative of an affiliation to the same group, whose correct name must be determined 
    ;in an additional step according to the specific rules of nomenclature that apply to the studied organisms.
    ;
    ;SeqID	MislabeledLevel	OriginalLabel	ProposedLabel	Confidence	OriginalTaxonomyPath	ProposedTaxonomyPath	PerRankConfidence
    END_SATIVA
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sativa: 0.9.3
    END_VERSIONS
    """
}

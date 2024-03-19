# ![nf-core/phylomarkercheck](docs/images/nf-core-phylomarkercheck_logo_light.png#gh-light-mode-only) ![nf-core/phylomarkercheck](docs/images/nf-core-phylomarkercheck_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/phylomarkercheck/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/phylomarkercheck/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/phylomarkercheck/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/phylomarkercheck/actions?query=workflow%3A%22nf-core+linting%22)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/phylomarkercheck/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/phylomarkercheck)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23phylomarkercheck-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/phylomarkercheck)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**biodiversitydata-se/sbdi-phylomarkercheck** is a bioinformatics pipeline that checks GTDB 16S sequences for phylogenetic signal with Sativa

Using Sativa [Kozlov et al. 2016], 16S sequences from GTDB are checked so that their phylogenetic signal is consistent with their taxonomy.

Before calling Sativa, sequences longer than 2000 nucleotides or containing Ns are removed, and the reverse complement of each is calculated.
Subsequently, sequences are aligned with HMMER [Eddy 2011] using the Barrnap [https://github.com/tseemann/barrnap] archaeal and bacterial 16S profiles respectively, 
and sequences containing more than 10% gaps are removed.
The remaining sequences are analyzed with Sativa, and sequences that are not phylogenetically consistent with their taxonomy are removed.

Files for the DADA2 (Callahan et al. 2016) methods `assignTaxonomy` and `addSpecies` are available, in three different versions each. 
The `assignTaxonomy` files contain taxonomy for domain, phylum, class, order, family, genus and species. 
(Note that it has been proposed that species assignment for short 16S sequences require 100% identity (Edgar 2018), so use species assignments with `assignTaxonomy` with caution.) 
The versions differ in the maximum number of genomes that we included per species: 1, 5 or 20, indicated by "n1", "n5" and "n20" in the file names respectively. 
Using the version with 20 genomes per species should increase the chances to identify an exactly matching sequence by the `addSpecies` algorithm, while using a file with many genomes 
per species could potentially give biases in the taxonomic annotations at higher levels by `assignTaxonomy`.
Our recommendation is hence to use the "n1" files for `assignTaxonomy` and "n20" for `addSpecies`.

All files are gzipped fasta files with 16S sequences, the assignTaxonomy associated with taxonomy hierarchies from domain to species whereas the `addSpecies` file have sequence identities and species names.
There are also fasta files with the original GTDB sequence names, with "correct" in their names.

Taxonomical annotation of 16S amplicons using this data is available as an optional argument to the nf-core/ampliseq Nextflow workflow from version 2.1: `--dada_ref_taxonomy sbdi-gtdb` 
(https://nf-co.re/ampliseq; Straub et al. 2020).

In addition to the fasta files, the workflow estimates phylogenetic trees from the original GTDB trees. 
As not all species in GTDB will have correct 16S sequences, the GTDB trees are first subset to contain only species for which the species representative genome has a correct 16S sequence.
Subsequently, branch lengths for the tree are optimized based on the original alignment of 16S sequences using IQTREE [Nguyen et al. 2015] with a GTR+F+I+G4 model.

## Usage

Create a parameter file, e.g. `params.yml`, similar to this:

```yml
markername: 'arc-ssu-r214'
input: 'input/arc-ssu-r214.fna'
hmm: 'https://raw.githubusercontent.com/tseemann/barrnap/master/db/arc.hmm'
hmmkey: '16S_rRNA'
gtdb_metadata: 'https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_metadata_r214.tsv.gz'
n_per_species: '1,5,20'
outdir: 'r214'
max_cpus: 12
non_gap_prop: 0.8 
phylogeny: 'https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_r214.tree'
model: 'GTR+F+I+G4'
```

And run the workflow like this:

```bash
nextflow run biodiversitydata-se/sbdi-phylomarkercheck \
   -profile <docker/singularity/.../institute> \
   --outdir <OUTDIR>
   -params-file params.yml
```

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/phylomarkercheck/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/phylomarkercheck/output).

## Credits

nf-core/phylomarkercheck was originally written by Daniel Lundin.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#phylomarkercheck` channel](https://nfcore.slack.com/channels/phylomarkercheck) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/phylomarkercheck for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

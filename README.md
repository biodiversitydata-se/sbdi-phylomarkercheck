# sbdi-phylomarkercheck

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

A set of directories under `<OUTDIR>` that you specified as argument for `--outdir` will be created.
Two of these are particularly interesting:
* `<OUTDIR>/correct`: Fasta files with 16S sequences with a taxonomically consistent phylogenetic signal
  - `<PREFIX>-n<N>.assignTaxonomy.fna.gz`: N sequences per species formatted for DADA2's `assignTaxonomy()`
  - `<PREFIX>-n<N>.addSpecies.fna.gz`: N sequences per species formatted for DADA2's `addSpecies()`
  - `<PREFIX>-n<N>.correct.fna.gz`: N sequences per species in GTDB's original format
* `<OUTDIR>/iqtree`: Tree files and taxonomy
  - `<PREFIX>-sprep.alnfna`: Fasta file with aligned 16S sequences
  - `<PREFIX>-sprep.brlenopt.treefile`: Newick formatted tree file
  - `<PREFIX>-sprep.taxonomy.tsv`: Taxonomy file for phylogenetic placement purposes
  - `<PREFIX>-sprep.brlenopt.iqtree`: IQTREE info file
  - `<PREFIX>-sprep.brlenopt.log`: Log file

## Credits

sbdi-phylomarkercheck was originally written by Daniel Lundin.

An earlier, manual, procedure to produce the corresponding files for releases up to r207 was designed in collaboration with Anders Andersson.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
 
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

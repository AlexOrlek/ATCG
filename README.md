ATCG (Alignment-based Tool for Comparative Genomics) is a command-line tool for pairwise comparison of nucleotide sequences; it is intended for small-to-medium-sized datasets (e.g. plasmids; bacterial genomes). [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) alignments are analysed to calculate various overall measures of sequence relatedness, which can then be visualised as dendrograms. The approach of ATCG is similar to that of a web-based tool called the genome-genome distance calculator ([GGDC](https://ggdc.dsmz.de/ggdc.php#)).

# Introduction

As input, ATCG can take sequence assemblies in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format; the assemblies can be any of the following:

* Complete genomes (e.g. plasmids; small collections of bacterial genomes)
* Complete metagenomes (e.g. sets of bacterial plasmids, where each set comprises plasmids from a single bacterial isolate)
* Incomplete genome/metagenome sequences (or a combination of complete/incomplete sequences). Incomplete sequences (contigs or scaffolds) must be affiliated to a known genome or metagenome and this affiliation is indicated in the sequence's FASTA header (see [Quick start](#quick-start))


# Requirements

* Linux
* [Python](https://www.python.org/) 2.7
* [SeqKit](https://github.com/shenwei356/seqkit)
* [bioawk](https://github.com/lh3/bioawk)
* [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (`blastn`)
* [R](https://www.r-project.org/) 3.3.1 or later with the following packages installed:
    * Bioconductor 3.4 and the [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) package
    * [gsubfn](https://cran.r-project.org/web/packages/gsubfn/index.html)


# Installation

```bash
git clone https://github.com/AlexOrlek/ATCG.git
```
# Quick start

The simplest way to use the tool is to provide a single multi-FASTA file with FASTA headers in the format `unit of analysis|subunit`. The `subunit` is only necessary where there are multiple subdivisions (multiple contigs or multiple genomes of a given metagenome).

If comparing genomes, FASTA headers might be as follows:<br>
genome1|contig1<br>
genome1|contig2<br>
genome2

If comparing plasmids by isolate, FASTA headers might be as follows:<br>
isolate1|plasmid1<br>
isolate1|plasmid2<br>
isolate2|plasmid1.contig1<br>
isolate2|plasmid1.contig2<br>
isolate2|plasmid2


The tool can be run as follows:

`python runpipeline.py -q genomes.fasta -o output-directory`


# Options and usage

Run `python runpipeline.py --help` to view a summary of all the options


# Output files

The below table shows the most important outputs

File/Directory                 | Description                                                                                       
------------------------------ | -------------------------------------------------------------------------------------------------
splitfastas/                   | directory containing FASTA files derived from the input multi-FASTA, split by unit of analysis i.e. genome (or metagenome)                                       
blast/			       | directory containing tsv files of blast alignments for each genome
output/distancestats.tsv       | columns of distance stats for each unique pairwise combination of genomes
output/dendrogram_[score].pdf  | dendrogram generated using a specified distance score column from the distancestats.tsv file
included.txt		       | names of genomes with detected blast alignments, that will therefore appear on the dendrogram
excluded.txt		       | names of any genomes with no detected blast alignnments (this file may well be blank)


# Background and methods

A methods paper will be written shortly. A brief outline is given below, and further information about the general approach can be found in a [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-60) by Meier-Kolthoff et al. describing the genome-genome distance calculator tool.

1. All-vs-all BLAST is conducted.
2. Overlapping alignments are trimmed.
3. For trimmed alignments, distance metrics are calculated; different metrics reflect different distance concepts: resemblance and containment. 
4. Using pairwise distances, a dendrogram is generated.


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)
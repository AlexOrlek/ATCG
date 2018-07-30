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
    * [purrr](https://github.com/tidyverse/purrr)


# Installation

```bash
git clone https://github.com/AlexOrlek/ATCG.git
```
# Quick start

Information in FASTA headers must be delineated using vetical bar(s) "|" and the headers should be in the format `unit of analysis|subunit`. The `subunit` is only necessary when indicating affiliation of e.g. contigs with genomes, or of bacterial plasmids with bacterial isolates. Additional information can be included in the header by delineating with additional vertical bar(s)

If comparing genomes, FASTA headers might be as follows:<br>
genome1|contig1<br>
genome1|contig2<br>
genome2<br>
genome3|contig1|additional information

If comparing plasmids by isolate, FASTA headers might be as follows:<br>
isolate1|plasmid1<br>
isolate1|plasmid2<br>
isolate2|plasmid1.contig1<br>
isolate2|plasmid1.contig2<br>
isolate2|plasmid2


For all-vs-all comparison, the tool can be run by providing a single multi-FASTA file using the -s flag; distance scores will be recorded and a dendrogram will be generated.

`python runpipeline.py -s genomes.fasta -o output-directory -t 8`


Alternatively, the tool can be run by providing 2 input files using flags -s1 and -s2; pairwise comparisons will be conducted between but not amongst sequence(s) in each file; a dendrogram will not be generated.

`python runpipeline.py -s1 query.fasta -s2 genomes.fasta -o output-directory -t 8`



# Options and usage

Run `python runpipeline.py --help` to view a summary of all the options


# Output files

The below table shows the most important outputs from running the pipeline with the -s input flag. Similar outputs are produced using -s1 and -s2 input flags. 

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

1. BLAST is conducted on assembled nucleotide sequences.
2. Overlapping alignments are trimmed.
3. For trimmed alignments, distance metrics are calculated; different metrics reflect different distance concepts: [resemblance and containment](https://www.cs.princeton.edu/courses/archive/spring13/cos598C/broder97resemblance.pdf). 
4. Using pairwise distances, a dendrogram is generated.


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)
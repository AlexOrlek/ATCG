ATCG (Alignment-based Tool for Comparative Genomics) is a command-line tool for comparison of nucleotide sequences; it is intended for small-to-medium-sized datasets (e.g. plasmids; small collections of bacterial genomes). [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) alignments are analysed to calculate various overall measures of pairwise sequence relatedness (expressed as distance scores). Pairwise distance scores can then be used to build a [distance-based](https://en.wikipedia.org/wiki/Distance_matrices_in_phylogeny) phylogenetic tree. The approach of ATCG is similar to that of a web-based tool called the genome-genome distance calculator ([GGDC](https://ggdc.dsmz.de/ggdc.php#)).

# Introduction

As input, ATCG can take nucleotide sequence assemblies in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format; the assemblies can comprise any of the following:

* Complete genomes (e.g. plasmids; small collections of bacterial genomes)
* Complete metagenomes (e.g. sets of bacterial plasmids, where each set comprises plasmids from a single bacterial isolate)
* Incomplete genome/metagenome sequences (or a combination of complete/incomplete sequences). Incomplete sequences (contigs or scaffolds) must be affiliated to a known genome or metagenome and this affiliation is indicated in the sequence's FASTA header (see [Quick start](#quick-start))

ATCG is appropriate if you want to:
* Compare sequences in terms of overall relatedness metrics (genome-genome distances, percentage identity, coverage breadth)
* Compare sequences in terms of their structural similarity (using the breakpoint distance metric; see [Henz _et al_. 2004](https://www.ncbi.nlm.nih.gov/pubmed/15166018))

ATCG is __not__ appropriate if you want to:
* Analyse many large genomes - alignment-based approaches are too time-consuming for this (instead, use locus-based typing methods such as [MLST](https://pubmlst.org/general.shtml) OR alignment-free comparative genomic tools such as [FastANI](https://github.com/ParBLiSS/FastANI))
* Generate a whole-genome multiple alignment (instead, use [progressiveMauve](http://darlinglab.org/mauve/user-guide/progressivemauve.html) or similar tool)


# Requirements

* Linux
* [Python](https://www.python.org/) 2.7 or Python 3 (tested on Python 3.5.2)
* [SeqKit](https://github.com/shenwei356/seqkit)
* [bioawk](https://github.com/lh3/bioawk)
* [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (`blastn`)
* [R](https://www.r-project.org/) 3.3.1 or later with the following packages installed:
    * [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html); [gsubfn](https://cran.r-project.org/web/packages/gsubfn/index.html); [purrr](https://github.com/tidyverse/purrr); [foreach](https://cran.r-project.org/web/packages/foreach/index.html); [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html); [data.table](https://cran.r-project.org/web/packages/data.table/index.html); [ape](https://cran.r-project.org/web/packages/ape/index.html)<br>

Run the following code in R to install the required R packages:<br>
```bash
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

install.packages("devtools")
library(devtools)
devtools::install_github("ggrothendieck/gsubfn")

install.packages("purrr")
install.packages("foreach")
install.packages("doParallel")
install.packages("data.table")
install.packages("ape")
```

    


# Installation

```bash
git clone https://github.com/AlexOrlek/ATCG.git
cd ATCG
```
# Quick start

Information in FASTA headers must be delineated using vertical bar(s) "|" and the headers should be in the format `unit of analysis|subunit`. The `subunit` is only necessary when indicating affiliation of e.g. contigs with genomes, or of bacterial plasmids with bacterial isolates. Additional information can be included in the header by delineating with additional vertical bar(s)

If comparing genomes, FASTA headers could be formatted as follows:<br>
genome1|contig1<br>
genome1|contig2<br>
genome2<br>
genome3|contig1|additional information

If comparing plasmids by isolate, FASTA headers could be formatted as follows:<br>
isolate1|plasmid1<br>
isolate1|plasmid2<br>
isolate2|plasmid1.contig1<br>
isolate2|plasmid1.contig2<br>
isolate2|plasmid2


For all-vs-all comparison, the tool can be run by providing a single multi-FASTA file using the `-s` flag; pairwise distance scores will be recorded and a phylogenetic tree will be generated.

`python runpipeline.py -s genomes.fasta -o output-directory -t 8`


Alternatively, if all-vs-all comparison is not required, 2 input (multi-)FASTA files can be provided using flags `-s1` and `-s2`; pairwise comparisons will be conducted between but not amongst sequence(s) in each file; the analysis will run faster than an all-vs-all comparison. Pairwise distances will be recorded but a tree will not be generated since distances are not available for all pairwise combinations.

`python runpipeline.py -s1 query.fasta -s2 genomes.fasta -o output-directory -t 8`



# Options and usage

Run `python runpipeline.py --help` to view a summary of all the options

By default, breakpoint distance is not calculated, but can be specified using the `--breakpoint` flag
By default, bootstrapping will not be conducted and the tree will therefore not show bootstrap confidence values. To conduct bootstrapping, the number of replicates is specified using the `-b` flag.

`python runpipeline.py -s genomes.fasta -o output-directory -t 8 -b 100 --breakpoint`

When bootstrapping is conducted, trimmed alignments will be resampled with replacement to produce replicate distance scores, from which replicate trees are produced, allowing confidence values to be shown on the original tree.


# Output files

The below table shows the most important outputs from running the pipeline with the -s input flag. Similar outputs are produced using -s1 and -s2 input flags, but trees will not be produced.

File/Directory         | Description                                                                                       
---------------------- | -------------------------------------------------------------------------------------------------
splitfastas/           | directory containing FASTA files derived from the input multi-FASTA, split by unit of analysis i.e. genome (or metagenome)                                       
blast/		       | directory containing tsv files of blast alignments for each genome
included.txt           | names of genomes with detected blast alignments, that will therefore appear in the distancestats.tsv file
excluded.txt	       | names of any genomes with no detected blast alignments (this file may well be blank)
output/		       | directory containing output files described below
distancestats.tsv      | columns of distance statistics for each unique pairwise combination of genomes
tree_[score].pdf       | tree generated using a specified distance score column from the distancestats.tsv file; plotted as a pdf
tree_[score].rds       | as above, but stored as an [rds file](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) which can be read, and the tree replotted 


If bootstrapping is specified, the following additional files will also be generated in the output directory:<br>
A file containing distance statistics for each bootstrap replicate (distancestats_bootstrapped.tsv)<br>
A pdf showing the original tree with bootstrap confidence values (tree\_[score]\_bootstrapped.pdf) is produced instead of tree_[score].pdf<br>
A list of replicate trees, that were used to calculate confidence values for the original tree, is stored as an rds file (tree\_[score]\_bootstrapped.rds)


# Background and methods

A paper describing the methods will be written shortly. A brief outline is given below, and further information about the general approach can be found in a [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-60) by Meier-Kolthoff _et al_. describing the genome-genome distance calculator tool.

1. BLAST is conducted on assembled nucleotide sequences.
2. Overlapping alignments are trimmed.
3. For trimmed alignments, distance metrics are calculated; different metrics reflect different distance concepts: [resemblance and containment](https://www.cs.princeton.edu/courses/archive/spring13/cos598C/broder97resemblance.pdf). Breakpoint distance can optionally be calculated.
4. If all-vs-all BLAST was run, then a tree is generated using a specified pairwise distance metric; optionally, the tree can be annotated with bootstrap confidence values, which are calculated by resampling trimmed alignments.


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)
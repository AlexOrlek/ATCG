ATCG (Alignment-based Tool for Comparative Genomics) is a command-line tool for comparison of nucleotide sequences; it is intended for small-to-medium-sized datasets (e.g. plasmids; small collections of bacterial genomes). [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) alignments are analysed to calculate various overall measures of pairwise sequence relatedness. The approach of ATCG is similar to that of a web-based tool called the genome-genome distance calculator ([GGDC](https://ggdc.dsmz.de/ggdc.php#)). However, ATCG provides more flexibility and additional analysis options such as assessment of structural similarity.

# Introduction

As input, ATCG can take nucleotide sequence assemblies in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format; the assemblies can comprise any of the following:

* Complete genomes (e.g. plasmids; small collections of bacterial genomes)
* Complete metagenomes (e.g. sets of bacterial plasmids, where each set comprises plasmids from a single bacterial isolate)
* Incomplete genome/metagenome sequences (or a combination of complete/incomplete sequences). Incomplete sequences (contigs or scaffolds) must be affiliated to a known genome or metagenome and this affiliation is indicated in the sequence's FASTA header (see [Input](#input))

ATCG is appropriate if you want to:
* Compare sequences in terms of overall relatedness metrics (genome-genome distances, percentage identity, coverage breadth)
* Use pairwise genome-genome distance scores to build a [distance-based](https://en.wikipedia.org/wiki/Distance_matrices_in_phylogeny) tree.
* Compare sequences in terms of their structural similarity (using the breakpoint distance metric; see [Henz _et al_. 2004](https://www.ncbi.nlm.nih.gov/pubmed/15166018))

ATCG is __not__ appropriate if you want to:
* Analyse many large genomes - alignment-based approaches are too time-consuming for this (instead, use locus-based typing methods such as [MLST](https://pubmlst.org/general.shtml) OR alignment-free comparative genomic tools such as [FastANI](https://github.com/ParBLiSS/FastANI))
* Generate a whole-genome multiple alignment (instead, use [progressiveMauve](http://darlinglab.org/mauve/user-guide/progressivemauve.html) or similar tool)


# Requirements

* Linux or MacOS (with the [Bash shell](https://en.wikibooks.org/wiki/Bash_Shell_Scripting#What_is_Bash?), which is the default shell on MacOS and many Linux distributions; tested using Bash versions 3.2, 4.1, 4.3)
* [Python](https://www.python.org/) 2.7 or Python 3 (tested using Python 3.5)
* [SeqKit](https://github.com/shenwei356/seqkit)
* [bioawk](https://github.com/lh3/bioawk)
* [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (`blastn`)
* [R](https://www.r-project.org/) 3.3.1 or later with the following packages installed:
    * [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html); [gsubfn](https://cran.r-project.org/web/packages/gsubfn/index.html); [purrr](https://github.com/tidyverse/purrr); [foreach](https://cran.r-project.org/web/packages/foreach/index.html); [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html); [data.table](https://cran.r-project.org/web/packages/data.table/index.html); [ape](https://cran.r-project.org/web/packages/ape/index.html)<br>

Run the following code in R to install the required R packages:<br>
```bash
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

install.packages("devtools",repo='https://cloud.r-project.org/')
library(devtools)
devtools::install_github("ggrothendieck/gsubfn")

install.packages("purrr",repo='https://cloud.r-project.org/')
install.packages("foreach",repo='https://cloud.r-project.org/')
install.packages("doParallel",repo='https://cloud.r-project.org/')
install.packages("data.table",repo='https://cloud.r-project.org/')
install.packages("ape",repo='https://cloud.r-project.org/')
```

    


# Installation

```bash
git clone https://github.com/AlexOrlek/ATCG.git
cd ATCG
```
# Input

Sequences are provided in FASTA format, a flexibly-defined format comprising header lines (prefixed by ">") and sequences. Here, we follow the common convention where the header line is permitted to have two parts, separated by a space: the identifier and an optional comment after the first space. Information in the header identifier must be delineated using vertical bar(s) "|" and adhere to the format: `unit of analysis|subunit`. The `subunit` is only necessary when indicating affiliation of e.g. contigs with genomes, or of bacterial plasmids with bacterial isolates.

If comparing genomes, FASTA headers could be formatted as follows:<br>
genome1|contig1<br>
genome1|contig2<br>
genome2<br>
genome3|contig1 additional information provided after the first space

If comparing isolates in terms of their overall plasmid genetic content, FASTA headers could be formatted as follows:<br>
isolate1|plasmid1<br>
isolate1|plasmid2<br>
isolate2|plasmid1<br>
isolate2|plasmid2

# Quick start

For all-vs-all comparison, the tool can be run by providing a single multi-FASTA file using the `-s` flag; pairwise distance scores will be recorded and a tree will be generated.

`python runpipeline.py -s genomes.fasta -o output-directory -t 8`


Alternatively, if all-vs-all comparison is not required, 2 input (multi-)FASTA files can be provided using flags `-s1` and `-s2`; pairwise comparisons will be conducted between but not amongst sequence(s) in each file; pairwise distances will be recorded but a tree will not be generated since distances are not available for all pairwise combinations. The order in which the FASTA files are provided to the -s1/-s2 flags will not affect results, however the sequences in the two files must be non-overlapping (ATCG will check this based on examination of the FASTA header identifiers).

`python runpipeline.py -s1 query.fasta -s2 genomes.fasta -o output-directory -t 8`



# Options and usage

`python runpipeline.py --help` produces a summary of all the options.

By default, the number of threads is 1, but multi-threading is recommended to reduce computing time; the number of threads to use is specified using the `-t` flag; the value must not exceed the number of threads available on your machine. 
By default, breakpoint distance and alignment length distribution statistics are not calculated; bootstrap confidence values are also not calculated by default.
To conduct bootstrapping, the number of replicates is specified using the `-b` flag; trimmed alignments will be resampled with replacement to produce replicate distance scores, from which replicate trees are produced, allowing bootstrap confidence values to be shown on the original tree.
Calculation of breakpoint distance (a measure of structural similarity) is specified using the `--breakpoint` flag.
Calculation of alignment length distribution statistics is specified using the `--alnlenstats` flag. The alignment length statistics provide information on the distribution of BLAST alignment lengths and are analogous to the widely used [assembly contiguity statistics](https://www.molecularecologist.com/2017/03/whats-n50/) e.g N50/L50.

`python runpipeline.py -s genomes.fasta -o output-directory -t 8 -b 100 --breakpoint --alnlenstats` runs the pipeline using 8 threads, with 100 bootstrap replicates, and calculation of breakpoint distance and alignment length distribution statistics.

To better understand the calculation of the various statistics, it may help to take a look at the [example](#Example) below.




# Output files

The below table shows the most important outputs from running the pipeline with the -s input flag. Similar outputs are produced using -s1 and -s2 input flags, but trees will not be produced.

File/Directory         | Description                                                                                       
---------------------- | -------------------------------------------------------------------------------------------------
splitfastas/           | directory containing FASTA files (and corresponding BLAST databases), derived from the input multi-FASTA, split by unit of analysis i.e. genome (or metagenome)
blast/		       | directory containing tsv files of blast alignments for each genome
included.txt           | names of genomes with detected blast alignments, that will therefore appear in the distancestats.tsv file
excluded.txt	       | names of any genomes with no detected blast alignments (this file may well be blank)
fastafilepaths.tsv     | genomes and corresponding FASTA file paths
blastdbfilepaths.tsv   | genomes and corresponding BLAST database file paths
seqlengths.tsv         | genomes and their lengths in bp
output/		       | directory containing output files described below
distancestats.tsv      | columns of distance statistics for each unique pairwise combination of genomes
tree_[score].pdf       | tree generated using a specified distance score column from the distancestats.tsv file; plotted as a pdf
tree_[score].rds       | as above, but stored as an [rds file](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) which can be read, and the tree replotted 
distobject_[score].rds | a "dist" object distance matrix derived from distancestats.tsv, stored as an rds file

If bootstrapping is specified, a pdf showing the original tree with bootstrap confidence values (tree\_[score]\_bootstrapped.pdf) is produced instead of tree_[score].pdf<br>
Also, the following additional files will be generated in the output directory:<br>
A file containing distance statistics for each bootstrap replicate (distancestats_bootstrapped.tsv)<br>
A list of replicate trees, that were used to calculate confidence values for the original tree, stored as an rds file (tree\_[score]\_bootstrapped.rds)


# Background and methods

A paper describing the methods will be written shortly, and further information about the general approach can be found in a [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-60) by Meier-Kolthoff _et al_. describing the similar genome-genome distance calculator tool ([GGDC](https://ggdc.dsmz.de/ggdc.php#)). A brief outline of the steps of ATCG is given below:

1. BLAST is conducted on assembled nucleotide sequences in both directions between each pair of genomes (i.e. genome A vs genome B and genome B vs genome A; that is, with genome A as the [query sequence](https://www.ncbi.nlm.nih.gov/books/NBK1734/) and genome B as the subject database sequence, and vice-versa).
2. Where alignment ranges overlap at the same region on the query genome or the subject genome, the shorter overlapping alignment is trimmed to eliminate the overlap. Trimming performed on the query/subject genome is applied to the corresponding range on the subject/query genome, accounting for the strand of the alignment. So, if an alignment is '-' strand (a reverse complement alignment), and the alignment range on the query sequence is trimmed from the 5' end, then the corresponding alignment range on the subject sequence will be trimmed by the same length on the 3' end. See the [example](#Example) below for a visual depiction of how alignment trimming works.
3. Using trimmed alignments from both BLAST directions, distance metrics are calculated; breakpoint distances and alignment length statistics can optionally be calculated. A description of the statistics produced by ATCG, and formulae for their calculation are given [here](misc/statistics_calculation.pdf).
4. If all-vs-all BLAST was run, then a tree is generated using a specified pairwise distance metric; optionally, the tree can be annotated with bootstrap confidence values, which are calculated by resampling trimmed alignments.


# Example

To clarify the calculation of statistics, the alignment trimming and statistics calculation stages of the ATCG pipeline can be run on very simple example BLAST output files, using the following code, called within the example folder:

`Rscript granges_example.R`

Results are produced in the example/output folder. Diagrams below show the alignments before trimming and after trimming; calculation of the statistics in the output folder can be done manually for the benefit of understanding, as is shown below for the calculation of percent identity. 

<p align="center">**Untrimmed alignments**</p>
<p align="center"><img src="example/images/untrimmed.JPG" alt="untrimmed" width="600"></p>
<p align="center">**Trimmed alignments**</p>
<p align="center"><img src="example/images/trimmed.JPG" alt="trimmed" width="600"></p>
<p align="center">**Calculation of percent identity from trimmed alignments**</p>
<p align="center"><img src="example/images/percent_identity.JPG" alt="percent identity calculation" width="600"></p>

Things to note:<br>
Notice how the red alignment is involved in alignment trimming: it overlaps with the turquoise alignment on sequence B and since the red alignment is longer, the turquoise alignment is trimmed. The red alignment also overlaps with the yellow alignment on sequence A; in this case, it is the shorter of the two alignments, so it is trimmed. However, because it is a '-' strand alignment, it is trimmed from the 3' end on sequence B.

The breakpoint distance is 1 meaning there were no pairs of sequences found to be adjacent and in the same relative order in sequence A and sequence B. While the red and yellow alignments are adjacent in both sequence A and B, one is '+' strand and the other is '-' strand so they are not in the same relative order.  

For the percent identity calculation, the numerator and denominator are multiplied by 2. This is because the numer of identical basepairs (numerator) and the alignment length (denominator) is the same in both BLAST directions. However, BLAST is not always symmetric; this is the rationale for aggregating across both BLAST directions. 

# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)
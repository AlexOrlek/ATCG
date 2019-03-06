#!/usr/bin/env python
import argparse, os, sys, subprocess
sourcedir=os.path.dirname(os.path.abspath(__file__))
cwdir=os.getcwd()
sys.path.append(sourcedir)

from pythonmods import runsubprocess

def positivefloat(x):
    x = float(x)
    if x < 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive value" %x)
    return x

def positiveint(x):
    x = int(x)
    if x < 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive value" %x)
    return x

def ltyint(x):
    x = int(x)
    if x < 0 or x > 6:
         raise argparse.ArgumentTypeError("%s is an invalid line type integer value" %x)
    return x


parser = argparse.ArgumentParser(description="ATCG: Alignment Based Tool for Comparative Genomics; visualise alignments",add_help=False)
#Help options
help_group = parser.add_argument_group('Help')
help_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')

#Input options                                                               
input_group = parser.add_argument_group('Input')
input_group.add_argument('-i','--inputdir', help='Input directory containing alignment files, and (optionally) tree file(s) (required)', required=True)
input_group.add_argument('-s','--syntax', help='Syntax of alignment file names and (optionally) annotation file names e.g. trimmedalignments_GENOMENAME.tsv GENOMENAME.gff (required)', nargs='+', required=True)
input_group.add_argument('-l','--sequencelengths', help='File specifying sequence lengths (required)', required=True)
input_group.add_argument('-c','--comparisons', help='File specifying which comparisons to visualise, output file name, and (optionally) annotation data/tree name (required)', required=True)
input_group.add_argument('-f','--features', help='Input directory containing feature annotation data files (optional)', required=False)

#Output options                                               
output_group = parser.add_argument_group('Output')
output_group.add_argument('-o','--out', help='Output directory (required)', required=True)

#Visualisation options                                  

#general
vis_group1 = parser.add_argument_group('Visualisation options: general')
vis_group1.add_argument('--comparisontype', help='Specifies whether comparisons are between a single reference and a set of queries, or if a chain of comparisons are to be visualised (default: chain)', default='chain', choices=['chain','singlereference'], type=str)
vis_group1.add_argument('--main', help='Title text (default: no title text)', default='NULL', type=str)
vis_group1.add_argument('--main_pos', help='Title text position if title text is provided (default: centre)', default='centre', choices=['centre','left','right'], type=str)
vis_group1.add_argument('--rightmargin', help='The margin to add to the right hand side, expressed as a proportion of the longest sequence length (default: 0.05)', default='0.05', type=positivefloat)
vis_group1.add_argument('--legend_orientation', help='Orientation of the legend (default: vertical)', default='vertical', choices=['vertical','horizontal'], type=str)
vis_group1.add_argument('--positivecols', help='Colours defining a colour scale reflecting percent identity of positive orientation alignments e.g. yellow red would plot low identity alignments yellow and high identity alignments red; more than 2 colours can be provided e.g. red orange yellow green blue violet (default: #f0cccc #b70000)', nargs='+',default=['#f0cccc', '#b70000'],type=str)
vis_group1.add_argument('--negativecols', help='As above but for negative orientation alignments (default: #ccccf0 #0000b7)', nargs='+',default=['#ccccf0', '#0000b7'],type=str)

#sequence-related
vis_group2 = parser.add_argument_group('Visualisation options: sequence-related')
vis_group2.add_argument('--sequencecols', help='Colour(s) of adjacent subsequences; default depends on --comparisontype; chain: "white"; singlereference: "light yellow" "light blue"', nargs='+',type=str)
vis_group2.add_argument('--dna_seg_labels', help='The names of the sequences; must provide same number of names as there are sequences. (default: use names of dna_segs if available)', default=['NULL'],nargs='+',type=str)
vis_group2.add_argument('--dna_seg_label_cex', help='The relative size of dna_seg_labels text (default: 0.8)', default='0.8',type=positivefloat)
vis_group2.add_argument('--dna_seg_label_col', help='The colour of dna_seg_labels text; provide either a single value or a space-separated list of the same length as the number of genomes (default: black)', default=['black'],nargs='+',type=str)
vis_group2.add_argument('--dna_seg_line', help='Should the line in the middle of the query segments be drawn; what colour(s) should be used; provide either true/false (case-insensitive) or colour; provide either a single value or a space-separated list of the same length as the number of genomes e.g. false would specify no line; blue would specify all lines as blue (default: false)', default=['false'],nargs='+',type=str)
vis_group2.add_argument('--minimum_gap_size', help='For chain reference comparison, the space between subsequences, expressed as the proportion of plotting width taken up by gaps (default: 0.05)', default='0.05',type=positivefloat)
vis_group2.add_argument('--scale', help='Should a global scale be shown at the bottom of the plot; provide true/false (case-insensitive) (default: false)', default='false',choices=['true','True','TRUE','false','False','FALSE'],type=str)
vis_group2.add_argument('--dna_seg_scale', help='Should per-genome scales be shown; provide true/false (case-insensitive) (default: true)', default='true',choices=['true','True','TRUE','false','False','FALSE'],type=str)
vis_group2.add_argument('--dna_seg_scale_cex', help='The relative size of per-genome scales (default: 0.7)', default='0.7',type=positivefloat)

#annotation-related
vis_group3 = parser.add_argument_group('Visualisation options: annotation-related')
vis_group3.add_argument('--annotation_gene_type', help='Specifies the default symbol to represent annotations (default: side_bars)', default='side_bars', choices=['arrowheads','arrows','headless_arrows','blocks','bars','points','text','lines','side_blocks','side_bars','side_points','side_text','side_lines','introns','exons','side_exons'],type=str)
vis_group3.add_argument('--annotation_outline_col', help='Specifies the default outline colour of the annotation symbol (default: light gray)', default='light gray',type=str)
vis_group3.add_argument('--annotation_fill_col', help='Specifies the default fill colour of the annotation symbol (default: black)', default='black',type=str)
vis_group3.add_argument('--annotation_lty', help='Specifies the default outline line type of the annotation symbol (default: 1)', default='1',type=ltyint)
vis_group3.add_argument('--annotation_lwd', help='Specifies the default outline line width of the annotation symbol (default: 0.8)', default='0.8',type=positivefloat)
vis_group3.add_argument('--annotationtxt_name', help='The tag name in the gff file to use for annotation text e.g. gene, product; unless specified, the gene tag will be used where available, and otherwise the product tag)', type=str)
vis_group3.add_argument('--annotationtxt_height', help='Specifies how much space to allocate to annotation text; calculated automatically based on annotation text length unless specified', type=positivefloat)
vis_group3.add_argument('--annotationtxt_rot', help='Specifies the annotation text rotation (default: 40)', default='40',type=int)
vis_group3.add_argument('--annotationtxt_cex', help='Specifies the relative size of annotation text (default: 0.5)', default='0.5',type=positivefloat)
vis_group3.add_argument('--annotationtxt_exclusion', help='Exclusion criteria for text annotations, command line argument(s) (optional)', nargs='+',type=str)
vis_group3.add_argument('--annotationtxt_inclusion', help='Inclusion criteria for text annotations, command line argument(s) (optional)', nargs='+',type=str)
vis_group3.add_argument('--annotationtxt_exclusion_file', help='Exclusion criteria for text annotations, text file (optional)')
vis_group3.add_argument('--annotationtxt_inclusion_file', help='Inclusion criteria for text annotations, text file (optional)')
vis_group3.add_argument('--annotationtxt_exclinc_casesensitive', action='store_true', help='If flag is provided, exclusion/inclusion criteria will be case-sensitive when matching annotations (default: case-insensitive)')
vis_group3.add_argument('--output_height', help='Specifies the height of the plot in inches; calculated automatically unless specified', type=positivefloat)
vis_group3.add_argument('--output_width', help='Specifies the width of the plot in inches; calculated automatically unless specified', type=positivefloat)


args = parser.parse_args()
outputpath=os.path.relpath(args.out, cwdir)

#!NB inclusionarg/exclusionsarg/syntaxvector must avoid comma

#check syntax flag correct
if 'GENOMENAME' not in str(args.syntax):
    print('Error: syntax flag must contain the placeholder text GENOMENAME which will be replaced with genome name specified in the comparison file')
    sys.exit()
if len(args.syntax)>2:
    print('Error: syntax flag must be provided with between 1 and 2 arguments; %i arguments provided'%len(args.syntax))
    sys.exit()

#sequencecol handling
if args.sequencecols==None:
    if args.comparisontype=='chain':
        args.sequencecols=["white"]
    else:
        args.sequencecols=["light yellow","light blue"]

#handling dna_seg_line - set true/false values to lowercase (will need to convert to TRUE/FALSE in R)
for i in range(len(args.dna_seg_line)):
    seglinearg=args.dna_seg_line[i]
    if seglinearg.lower()=='true' or seglinearg.lower()=='false':
        args.dna_seg_line[i]=seglinearg.lower()

#handling flags with auto setting
if args.annotationtxt_name==None:
    args.annotationtxt_name='auto'
if args.annotationtxt_height==None:
    args.annotationtxt_height='auto'
if args.output_height==None:
    args.output_height='auto'
if args.output_width==None:
    args.output_width='auto'
if args.output_width==None:
    args.output_width='auto'


#annotation text handling - only allow either command line or file for inclusion/exclusion criteria
if args.annotationtxt_exclusion!=None and args.annotationtxt_exclusion_file!=None:
    print('Error: you can provide either --annotatointxt_exclusion or --annotatointxt_exclusion_file arguments, but not both')
    sys.exit()
if args.annotationtxt_inclusion!=None and args.annotationtxt_inclusion_file!=None:
    print('Error: you can provide either --annotatointxt_exclusion or --annotatointxt_exclusion_file arguments, but not both')
    sys.exit()

if args.annotationtxt_exclusion==None and args.annotationtxt_exclusion_file==None:
    exclusionpresent='exclusionabsent'
    exclusionarg='placeholder'
elif args.annotationtxt_exclusion!=None:
    exclusionpresent='commandline'
    exclusionarg=str(','.join(args.annotationtxt_exclusion))
else:
    exclusionpresent='filepath'
    exclusionarg=args.annotationtxt_exclusion_file
    if os.path.isfile(exclusionarg)==False:
        print('Error: %s is not a valid filepath'%exclusionarg)
        sys.exit()

if args.annotationtxt_inclusion==None and args.annotationtxt_inclusion_file==None:
    inclusionpresent='inclusionabsent'
    inclusionarg='placeholder'
elif args.annotationtxt_inclusion!=None:
    inclusionpresent='commandline'
    inclusionarg=str(','.join(args.annotationtxt_inclusion))
else:
    inclusionpresent='filepath'
    inclusionarg=args.annotationtxt_inclusion_file
    if os.path.isfile(inclusionarg)==False:
        print('Error: %s is not a valid filepath'%inclusionarg)
        sys.exit()


#handle filepaths to directory
args.inputdir=str(args.inputdir).rstrip('/')


runsubprocess(['mkdir -p %s'%outputpath],shell=True)

if args.features==None:
    runsubprocess(['Rscript','%s/genoplotr.R'%sourcedir,str(args.inputdir),str(','.join(args.syntax)),str(args.sequencelengths),str(args.comparisons),outputpath,str(args.comparisontype),str(args.main),str(args.main_pos),str(args.rightmargin),str(','.join(args.sequencecols)),str(','.join(args.dna_seg_labels)),str(args.dna_seg_label_cex),str(','.join(args.dna_seg_label_col)),str(','.join(args.dna_seg_line)),str(args.minimum_gap_size),str(args.scale).lower(),str(args.dna_seg_scale).lower(),str(args.dna_seg_scale_cex),str(args.output_height),str(args.output_width),str(args.legend_orientation),','.join(args.positivecols),','.join(args.negativecols),sourcedir,'featuresabsent'])
else:
    runsubprocess(['Rscript','%s/genoplotr.R'%sourcedir,str(args.inputdir),str(','.join(args.syntax)),str(args.sequencelengths),str(args.comparisons),outputpath,str(args.comparisontype),str(args.main),str(args.main_pos),str(args.rightmargin),str(','.join(args.sequencecols)),str(','.join(args.dna_seg_labels)),str(args.dna_seg_label_cex),str(','.join(args.dna_seg_label_col)),str(','.join(args.dna_seg_line)),str(args.minimum_gap_size),str(args.scale).lower(),str(args.dna_seg_scale).lower(),str(args.dna_seg_scale_cex),str(args.output_height),str(args.output_width),str(args.legend_orientation),','.join(args.positivecols),','.join(args.negativecols),sourcedir,'featurespresent',str(args.features).rstrip('/'),str(args.annotationtxt_name),str(args.annotationtxt_height),str(args.annotationtxt_rot),str(args.annotationtxt_rot),exclusionpresent,exclusionarg,inclusionpresent,inclusionarg,str(args.annotationtxt_exclinc_casesensitive),str(args.annotation_gene_type),str(args.annotation_outline_col),str(args.annotation_fill_col),str(args.annotation_lty),str(args.annotation_lwd)])
print('finished plotting visualisation')



#OLD
              
#input_group.add_argument('--tree', help='Tree file, used to plot tree alongside comparisons (optional)')
#vis_group.add_argument('--annotation_cex', help='Specifies the default cex (scaling) of the annotation symbol (default: 1)', default='1',type=positiveint)
#vis_group.add_argument('--annotation_cex', help='Specifies the default cex (scaling) of the annotation symbols; calculated automatically based on output dimensions unless specified, or provided on a per-gene basis in the gff file (optional)', type=positivefloat)
#vis_group.add_argument('--annotation_lwd', help='Specifies the default height of the annotation symbols; calculated automatically based on output dimensions unless specified, or provided on a per-gene basis in the gff file (optional)',type=positivefloat)
#if args.annotation_lwd==None:
#    args.annotation_lwd='auto'

#, donotplot; if donotplot, annotation text will not be plotted

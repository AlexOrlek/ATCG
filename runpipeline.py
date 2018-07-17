import argparse, os, sys, time
from mymod import runsubprocess


def proportion(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x
            

parser = argparse.ArgumentParser(description='Run pipeline scripts')
parser.add_argument('-o','--out', help='Output directory (required)', required=True)
parser.add_argument('-q','--plasmidqueries', help='Plasmid sequences for all samples, formatted as a multi-FASTA file (required)', required=True)
parser.add_argument('-d','--plasmiddbs', help='Text file containing filepaths to BLAST databases for each sample (required)', required=True)
parser.add_argument('-a','--analysislevel', help='Is the analysis conducted at the "individual" plasmid level or "samplelevel" (required)', required=True)
parser.add_argument('-e','--evalue', help='BLAST e-value cutoff (default: 1e-8)', default=1e-8, type=float) #1e-8 is used in ggdc web pipeline - see Meier-Kolthoff 2014
parser.add_argument('-w','--wordsize', help='BLAST word size (default: 38)', default=38, type=int) #38 is used in ggdc web pipleine?
parser.add_argument('-t','--threads', help='Number of threads to use (default: 8)', default=8, type=int)
parser.add_argument('-s','--distscore', help='Distance score to use (default: DistanceScore_d8d9)', default='DistanceScore_d8d9', type=str)
args = parser.parse_args()


outputpath=os.path.relpath(args.out, os.path.dirname(os.path.abspath(__file__)))

start=time.time()
runsubprocess(['python','runblast.py',outputpath, str(args.plasmidqueries), str(args.plasmiddbs), str(args.evalue), str(args.wordsize), str(args.threads)])
later=time.time()
print(later-start, 'runtime; finished running blast')
runsubprocess(['bash','splitblast.sh',outputpath,str(args.plasmiddbs)])
later=time.time()
print(later-start, 'runtime; finished splitting alignments')
runsubprocess(['bash','getseqlengths.sh',outputpath, str(args.plasmidqueries)])
later=time.time()
print(later-start, 'runtime; finished getting sequence lengths')
if str(args.analysislevel)=='individual':
    runsubprocess(['Rscript','granges.R',outputpath, str(args.threads)])
elif str(args.analysislevel)=='samplelevel':
    runsubprocess(['Rscript','granges_samplelevel.R',outputpath, str(args.threads)])
else:
    print('unrecognised samplelevel argument')
    sys.exit()
later=time.time()
print(later-start, 'runtime; finished trimming alignments')
if str(args.distscore)=='DistanceScore_d8d9':
    runsubprocess(['Rscript','dendrogram.R',outputpath, 'DistanceScore_d8'])
    later=time.time()
    print(later-start, 'runtime; finished plotting dendrogram')
    runsubprocess(['Rscript','dendrogram.R',outputpath, 'DistanceScore_d9'])
    later=time.time()
    print(later-start, 'runtime; finished plotting dendrogram')
else:
    runsubprocess(['Rscript','dendrogram.R',outputpath, str(args.distscore)])
    later=time.time()
    print(later-start, 'runtime; finished plotting dendrogram')

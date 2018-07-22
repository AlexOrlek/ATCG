import argparse, os, sys, time
sys.path.append('./')
from pythonmods import runsubprocess


def proportion(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x
            

parser = argparse.ArgumentParser(description='Run pipeline scripts')
parser.add_argument('-q','--queries', help='Sequences, formatted as a multi-FASTA file (required)', required=True)
parser.add_argument('-o','--out', help='Output directory (required)', required=True)
parser.add_argument('-f','--fastafiles', help='Text file containing filepaths to fasta files for each genome (default: fastafilepaths.tsv created in output directory)', required=False)
parser.add_argument('-d','--blastdbs', help='Text file containing filepaths to BLAST databases for each genome (default: blastdbfilepaths.tsv created in output directory)', required=False)
parser.add_argument('-e','--evalue', help='BLAST e-value cutoff (default: 1e-8)', default=1e-8, type=float) #1e-8 is used in ggdc web pipeline - see Meier-Kolthoff 2014
parser.add_argument('-w','--wordsize', help='BLAST word size (default: 38)', default=38, type=int) #38 is used in ggdc web pipleine?
parser.add_argument('-t','--threads', help='Number of threads to use (default: 8)', default=8, type=int)
parser.add_argument('-s','--distscore', help='Distance score to use to construct dendrogram (can specify multiple parameters; default: DistanceScore_d8 DistanceScore_d9)', nargs='+', default=['DistanceScore_d8', 'DistaneScore_d9'], choices=['DistanceScore_d0','DistanceScore_d4','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9'], type=str)
args = parser.parse_args()
outputpath=os.path.relpath(args.out, os.path.dirname(os.path.abspath(__file__)))



start=time.time()
if args.fastafiles==None and args.blastdbs==None:
    args.fastafiles='%s/fastafilepaths.tsv'%outputpath
    args.blastdbs='%s/blastdbfilepaths.tsv'%outputpath
    runsubprocess(['bash','splitfasta.sh',outputpath, str(args.queries), str(args.threads)])
    runsubprocess(['python','renamefastas.py',outputpath, str(args.queries)])
    runsubprocess(['python','makeblastdbs.py',outputpath, str(args.fastafiles)])
    later=time.time()
    print(later-start, 'runtime; finished creating blast databases')
if args.fastafiles!=None and args.blastdbs==None:
    args.blastdbs='%s/blastdbfilepaths.tsv'%outputpath
    runsubprocess(['python','makeblastdbs.py',outputpath, str(args.fastafiles)])   
    later=time.time()
    print(later-start, 'runtime; finished creating blast databases')
runsubprocess(['python','runblast.py',outputpath, str(args.queries), str(args.blastdbs), str(args.evalue), str(args.wordsize), str(args.threads)])
later=time.time()
print(later-start, 'runtime; finished running blast')
runsubprocess(['bash','reformatblastoutput.sh',outputpath,str(args.blastdbs)])
later=time.time()
print(later-start, 'runtime; finished splitting alignments')
runsubprocess(['bash','getseqlengths.sh',outputpath, str(args.queries)])
later=time.time()
print(later-start, 'runtime; finished getting sequence lengths')
runsubprocess(['Rscript','granges.R',outputpath, str(args.threads)])
later=time.time()
print(later-start, 'runtime; finished trimming alignments')
dendargs=['Rscript','dendrogram.R',outputpath]
dendargs.extend(args.distscore)
runsubprocess(dendargs)
later=time.time()
print(later-start, 'runtime; finished plotting dendrogram using distance score(s): %s'%arg.distscore)
                    



#OLD
#parser.add_argument('-a','--analysislevel', help='Is the analysis conducted at the "individual" plasmid level or "samplelevel" (required)', required=True)

# if str(args.analysislevel)=='individual':
#     runsubprocess(['Rscript','granges.R',outputpath, str(args.threads)])
# elif str(args.analysislevel)=='samplelevel':
#else:
#    print('unrecognised samplelevel argument')
#    sys.exit()



# if str(args.distscore)=='DistanceScore_d8d9':
#     runsubprocess(['Rscript','dendrogram.R',outputpath, 'DistanceScore_d8'])
#     later=time.time()
#     print(later-start, 'runtime; finished plotting dendrogram')
#     runsubprocess(['Rscript','dendrogram.R',outputpath, 'DistanceScore_d9'])
#     later=time.time()
#     print(later-start, 'runtime; finished plotting dendrogram')
# else:
#     runsubprocess(['Rscript','dendrogram.R',outputpath, str(args.distscore)])
#     later=time.time()
#     print(later-start, 'runtime; finished plotting dendrogram')

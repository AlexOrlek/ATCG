#!/usr/bin/env python
import argparse, os, sys, time,subprocess, signal
sourcedir=os.path.dirname(os.path.abspath(__file__))
cwdir=os.getcwd()
sys.path.append(sourcedir)

from pythonmods import runsubprocess


def default_sigpipe():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def positiveint(x):
    x = int(x)
    if x <= 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" %x)
    return x

def runtime():
    try:
        runtime=time.perf_counter()
    except:
        runtime=time.time()
    return(float(runtime))


parser = argparse.ArgumentParser(description='Run pipeline scripts')
parser.add_argument('-s','--sequences', help='Sequences, for all-vs-all pairwise comparison', required=False)
parser.add_argument('-s1','--sequences1', help='First set of sequence(s), for pairwise comparison against second set', required=False)
parser.add_argument('-s2','--sequences2', help='Second set of sequence(s), for pairwise comparison against first set', required=False)
parser.add_argument('-o','--out', help='Output directory (required)', required=True)
parser.add_argument('-e','--evalue', help='BLAST e-value cutoff (default: 1e-8)', default=1e-8, type=float) #1e-8 is used in ggdc web pipeline - see Meier-Kolthoff 2014
parser.add_argument('-w','--wordsize', help='BLAST word size (default: 38)', default=38, type=int) #38 is used in ggdc web pipleine?
parser.add_argument('-t','--threads', help='Number of threads to use (default: 1)', default=1, type=int)
parser.add_argument('-b','--boot', help='Number of bootstraps to run (default: no bootstrapping)', default=0, type=positiveint)
parser.add_argument('-d','--distscore', help='Distance score to use to construct tree (can specify multiple parameters; default: DistanceScore_d8 DistanceScore_d9)', nargs='+', default=['DistanceScore_d8', 'DistanceScore_d9'], choices=['DistanceScore_d0','DistanceScore_d4','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9'], metavar='',type=str)
parser.add_argument('-m','--treemethod', help='Tree building method (can specify dendrogram and/or phylogeny, or none (i.e. no tree will be plotted); default: dendrogram)', nargs='+', default=['dendrogram'], choices=['dendrogram','phylogeny','none'], metavar='',type=str)
parser.add_argument('--breakpoint', action='store_true', help='Calculate breakpoint statistics (default: do not calculate)')
parser.add_argument('--alnlenstats', action='store_true', help='Calculate alignment length distribution statistics (default: do not calculate)')
args = parser.parse_args()
outputpath=os.path.relpath(args.out, cwdir)


startruntime=runtime()
if args.sequences==None and args.sequences1==None and args.sequences2==None:
    parser.error('as input, you must either provide --sequences or both --sequences1 and --sequences2')
if args.sequences!=None and args.sequences1!=None and args.sequences2!=None:
    parser.error('as input, you must either provide --sequences or both --sequences1 and --sequences2')
if (args.sequences1==None and args.sequences2!=None) or (args.sequences1!=None and args.sequences2==None):
    parser.error('as input, you must either provide --sequences or both --sequences1 and --sequences2')


if args.sequences!=None:
    blasttype='allvallpairwise'
    fastadir='%s/splitfastas'%outputpath
    fastafiles='%s/fastafilepaths.tsv'%outputpath
    blastdbs='%s/blastdbfilepaths.tsv'%outputpath
    runsubprocess(['bash','%s/splitfasta.sh'%sourcedir,outputpath, fastadir, str(args.sequences),str(args.threads)])
    runsubprocess(['python','%s/renamefastas.py'%sourcedir,outputpath, fastadir, str(args.sequences),fastafiles,blastdbs])
    runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,fastadir, fastafiles, str(args.threads),sourcedir])
    laterruntime=runtime()
    print(laterruntime-startruntime, 'runtime; finished creating blast databases')
    p=subprocess.Popen(['bash','%s/runblast.sh'%sourcedir,outputpath, str(args.sequences), blastdbs, str(args.evalue), str(args.wordsize), str(args.threads)], preexec_fn=default_sigpipe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr= p.communicate()
    try:
        print('{} {}'.format(stdout.decode(), 'stdout'))
    except:
        pass
    try:
        print('{} {}'.format(stderr.decode(), 'stderr'))
    except:
        pass
    if p.returncode!=0:
        sys.exit()
    laterruntime=runtime()
    print(laterruntime-startruntime, 'runtime; finished running blast')
    runsubprocess(['bash','%s/reformatblastoutput.sh'%sourcedir,outputpath,blastdbs,sourcedir])
    laterruntime=runtime()
    print(laterruntime-startruntime, 'runtime; finished reformatting alignments')
    runsubprocess(['bash','%s/getseqlengths.sh'%sourcedir,outputpath,blasttype,str(args.sequences)])
    laterruntime=runtime()
    print(laterruntime-startruntime, 'runtime; finished getting sequence lengths')
    
if args.sequences==None:
    blasttype='pairwise'
    fastadir1='%s/splitfastas1'%outputpath
    fastadir2='%s/splitfastas2'%outputpath
    fastafiles1='%s/fastafilepaths1.tsv'%outputpath
    blastdbs1='%s/blastdbfilepaths1.tsv'%outputpath
    fastafiles2='%s/fastafilepaths2.tsv'%outputpath
    blastdbs2='%s/blastdbfilepaths2.tsv'%outputpath
    runsubprocess(['bash','%s/splitfasta.sh'%sourcedir,outputpath, fastadir1, str(args.sequences1),str(args.threads)])
    runsubprocess(['bash','%s/splitfasta.sh'%sourcedir,outputpath, fastadir2, str(args.sequences2),str(args.threads)])
    runsubprocess(['python','%s/renamefastas.py'%sourcedir,outputpath, fastadir1, str(args.sequences1), fastafiles1,blastdbs1])
    runsubprocess(['python','%s/renamefastas.py'%sourcedir,outputpath, fastadir2, str(args.sequences2), fastafiles2,blastdbs2])
    #check there is no overlap between fastas provided in -s1 and -s2
    s1dir=os.listdir('%s/splitfastas1'%outputpath)
    s2dir=os.listdir('%s/splitfastas2'%outputpath)
    s1fastas=[f for f in s1dir if f.endswith('.fasta')]
    s2fastas=[f for f in s2dir if f.endswith('.fasta')]
    overlap=len(set(s1fastas).intersection(set(s2fastas)))
    if overlap>0:
        parser.error('there must be no overlap between fasta identifiers contained within the fasta files provided using the -s1 and -s2 flags')
    runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,fastadir1, fastafiles1, str(args.threads)])
    runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,fastadir2, fastafiles2, str(args.threads)])
    laterruntime=runtime()
    print(laterruntime-startruntime, 'runtime; finished creating blast databases')
    p=subprocess.Popen(['bash','%s/runblast.sh'%sourcedir,outputpath, str(args.sequences1), blastdbs2, str(args.evalue), str(args.wordsize), str(args.threads)], preexec_fn=default_sigpipe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr= p.communicate()
    try:
        print('{} {}'.format(stdout.decode(), 'stdout'))
    except:
        pass
    try:
        print('{} {}'.format(stderr.decode(), 'stderr'))
    except:
        pass
    if p.returncode!=0:
        sys.exit()
    p=subprocess.Popen(['bash','runblast.sh',outputpath, str(args.sequences2), blastdbs1, str(args.evalue), str(args.wordsize), str(args.threads)], preexec_fn=default_sigpipe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr= p.communicate()
    try:
        print('{} {}'.format(stdout.decode(), 'stdout'))
    except:
        pass
    try:
        print('{} {}'.format(stderr.decode(), 'stderr'))
    except:
        pass
    if p.returncode!=0:
        sys.exit()
        
    laterruntime=runtime()
    print(laterruntime-startruntime, 'runtime; finished running blast')
    
    blastdbs='%s/blastdbfilepaths_combined.tsv'%outputpath
    runsubprocess(['cat %s %s | sort -k1,1V > %s'%(blastdbs1,blastdbs2,blastdbs)],shell=True)

    runsubprocess(['bash','%s/reformatblastoutput.sh'%sourcedir,outputpath,blastdbs])
    laterruntime=runtime()
    print(laterruntime-startruntime, 'runtime; finished reformatting alignments')
    runsubprocess(['bash','%s/getseqlengths.sh'%sourcedir,outputpath,blasttype,str(args.sequences1),str(args.sequences2)])
    laterruntime=runtime()
    print(laterruntime-startruntime, 'runtime; finished getting sequence lengths')

runsubprocess(['Rscript','%s/granges.R'%sourcedir,outputpath, str(args.threads), str(args.breakpoint),str(args.alnlenstats),str(args.boot)])
laterruntime=runtime()
print(laterruntime-startruntime, 'runtime; finished trimming alignments')
if blasttype=='allvallpairwise':
    if 'none' not in args.treemethod:
      if 'dendrogram' in args.treemethod:
          dendrogram='True'
      else:
          dendrogram='False'
      if 'phylogeny' in args.treemethod:
          phylogeny='True'
      else:
          phylogeny='False'
          
      linecount=int(runsubprocess(['wc -l < %s/included.txt'%outputpath],shell=True))
      if linecount < 3:
          print('Warning: there are only %i samples: a tree cannot be constructed'%linecount)
      else:
          treeargs=['Rscript','%s/tree.R'%sourcedir,outputpath,str(args.threads),str(args.boot),dendrogram,phylogeny]
          treeargs.extend(args.distscore)
          runsubprocess(treeargs)
          laterruntime=runtime()
          print(laterruntime-startruntime, 'runtime; finished plotting tree(s) using distance score(s): %s'%args.distscore)
#!/usr/bin/env python
import argparse, os, sys, time,subprocess, signal
from argparse import RawTextHelpFormatter

sourcedir=os.path.dirname(os.path.abspath(__file__))
cwdir=os.getcwd()
sys.path.append(sourcedir)

from pythonmods import runsubprocess,splitfastas


def default_sigpipe():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def positiveint(x):
    x = int(x)
    if x < 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" %x)
    return x

def runtime():
    try:
        runtime=time.perf_counter()
    except:
        runtime=time.time()
    return(float(runtime))


parser = argparse.ArgumentParser(description="ATCG: Alignment Based Tool for Comparative Genomics",add_help=False,formatter_class=RawTextHelpFormatter)
#Help options
help_group = parser.add_argument_group('Help')
help_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
#Input options                                                               
input_group = parser.add_argument_group('Input')
input_group.add_argument('-s','--sequences', help='Sequences, for all-vs-all pairwise comparison (required if -1 and -2 flags not provided)', required=False,nargs='?',type=argparse.FileType('r'),default=sys.stdin)
input_group.add_argument('-1','--sequences1', help='First set of sequence(s), for pairwise comparison against second set (required if -s flag not provided)', required=False,nargs='?',type=argparse.FileType('r'),default=sys.stdin)
input_group.add_argument('-2','--sequences2', help='Second set of sequence(s), for pairwise comparison against first set (required if -s flag not provided)', required=False,nargs='?',type=argparse.FileType('r'),default=sys.stdin)
#Output options                                               
output_group = parser.add_argument_group('Output')
output_group.add_argument('-o','--out', help='Output directory (required)', required=True)
output_group.add_argument('-d','--distscore', help='Distance score to use to construct tree; can specify multiple parameters (default: DistanceScore_d8 DistanceScore_d9)', nargs='+', default=['DistanceScore_d8', 'DistanceScore_d9'], choices=['DistanceScore_d0','DistanceScore_d4','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9'], metavar='', type=str)
output_group.add_argument('-m','--treemethod', help='Tree building method; can specify dendrogram and/or phylogeny, or none (i.e. no tree will be plotted); (default: dendrogram)', nargs='+', default=['dendrogram'], choices=['dendrogram','phylogeny','none'], metavar='', type=str)
output_group.add_argument('--breakpoint', action='store_true', help='If flag is provided, calculate breakpoint statistics (default: do not calculate)')
output_group.add_argument('--alnlenstats', action='store_true', help='If flag is provided, calculate alignment length distribution statistics (default: do not calculate)')
output_group.add_argument('--trimmedalignments', action='store_true', help='If flag is provided, write to file trimmed alignments for each sample (default: do not output)')
output_group.add_argument('--disjointalignments', action='store_true', help='If flag is provided, write to file disjoint alignments for each sample (default: do not output)')
#BLAST options                                  
blast_group = parser.add_argument_group('BLAST options')
blast_group.add_argument('--evalue', help='BLAST e-value cutoff (default: 1e-8)', default=1e-8, type=float) #1e-8 is used in ggdc web pipeline - see Meier-Kolthoff 2014
blast_group.add_argument('--wordsize', help='BLAST word size (ATCG default for blastn: 38; ATCG default for dc-megablast: 12; ATCG default for megablast: 28)', type=int) #38 is used in ggdc web pipleine?                 
blast_group.add_argument('--task', help='BLAST task (default: blastn)', default='blastn', choices=['blastn','dc-megablast','megablast'], type=str)
blast_group.add_argument('--bidirectionalblast', action='store_true', help='If flag is provided, BLAST is conducted in both directions and results are averaged (slower runtime) (default: conduct BLAST in one direction only)')
#Alignment filtering options                                                   
alignment_group = parser.add_argument_group('Alignment filtering options')
alignment_group.add_argument('-l','--lengthfilter', help='Length threshold (in basepairs) used to filter trimmed alignments prior to calculating breakpoint distance and alignment length statistics (default: 100)', default=100, type=positiveint)
alignment_group.add_argument('-r','--alnrankmethod', help='Parameter used for selecting best alignment (default: bitscore)', default='bitscore', choices=['bitscore','alnlen','pid'], type=str)
#Other options
other_group = parser.add_argument_group('Other')
other_group.add_argument('-t','--threads', help='Number of threads to use (default: 1)', default=1, type=int)
other_group.add_argument('-b','--boot', help='Number of bootstraps to run (default: no bootstrapping)', default=0, type=positiveint)
output_group.add_argument('--blastonly', action='store_true', help='If flag is provided, only output pairwise blast results; do not calculate distance statistics (default: run full pipeline)')
output_group.add_argument('--keep', default=1, choices=[0,1,2], type=int, help='Level of file retention (default: 1)\n  '
                          ' 0 = do not retain splitfastas or blast directories,\n'
                          '   1 = do not retain splitfastas directory,\n'
                          '   2 = retain all intermediate output')

args = parser.parse_args()
outputpath=os.path.relpath(args.out, cwdir)

startruntime=runtime()

#check sequence input flags used correctly (if flags not provided, name of file handle will be <stdin> i.e. empty stdin)
if args.sequences.name=='<stdin>' and args.sequences1=='<stdin>' and args.sequences2=='<stdin>':
    parser.error('as input, you must either provide --sequences or both --sequences1 and --sequences2')
if args.sequences.name!='<stdin>' and args.sequences1.name!='<stdin>' and args.sequences2.name!='<stdin>':
    parser.error('as input, you must either provide --sequences or both --sequences1 and --sequences2')
if (args.sequences1.name=='<stdin>' and args.sequences2.name!='<stdin>') or (args.sequences1.name!='<stdin>' and args.sequences2=='<stdin>'):
    parser.error('as input, you must either provide --sequences or both --sequences1 and --sequences2')

#set default word sizes if none provided; check word size is within allowed limits for the blast task
if args.wordsize==None:
    if args.task=='blastn':
        args.wordsize=int(38)
    if args.task=='dc-megablastn':
        args.wordsize=int(12)
    if args.task=='megablast':
        args.wordsize=int(28)
else:
    if args.task=='dc-megablast':
        if args.wordsize>int(12) or args.wordsize<int(11):
            parser.error('if using dc-megablast, word size must be either 11 or 12')
    if args.task=='blastn':
        if args.wordsize<int(4):
            parser.error('if using blastn, word size must be at least 4, and a higher value would be sensible, to reduce computation time when running genome-genome searches')
    if args.task=='megablast':
        if args.wordsize<int(4):
            parser.error('if using megablast, word size must be at least 4, and a higher value would be sensible, to reduce computation time when running genome-genome searches')

#check --keep flag and --blastonly flag are not used in combination (makes no sense)
if args.keep==0 and args.blastonly==True:
    parser.error('combining --keep 0 and --blastonly options will result in no final output being produced')


          
noblasthits=False

if args.sequences.name!='<stdin>':
    blasttype='allvallpairwise'
    fastadir='%s/splitfastas'%outputpath
    fastafiles='%s/fastafilepaths.tsv'%outputpath
    blastdbs='%s/blastdbfilepaths.tsv'%outputpath
    if os.path.exists(fastadir):
        sys.exit('Error: %s directory already exists, delete directory and try again'%fastadir)
    splitfastas(args.sequences,fastadir,fastafiles,blastdbs)
    runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,fastafiles,str(args.threads),sourcedir])
    laterruntime=runtime()
    #print(laterruntime-startruntime, 'runtime; finished creating blast databases')
    print('finished creating blast databases')
    p=subprocess.Popen(['bash','%s/runblast.sh'%sourcedir,outputpath, fastadir, blastdbs, str(args.evalue), str(args.wordsize), str(args.task),str(args.threads),str(args.bidirectionalblast)], preexec_fn=default_sigpipe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr= p.communicate()
    #try:
    #    print('{} {}'.format(stdout.decode(), 'stdout'))
    #except:
    #    pass
    #try:
    #    print('{} {}'.format(stderr.decode(), 'stderr'))
    #except:
    #    pass
    if p.returncode!=0:
        sys.exit('unexpected error when running runblast.sh')
    laterruntime=runtime()
    #print(laterruntime-startruntime, 'runtime; finished running blast')
    print('finished running blast')
    p=subprocess.Popen(['bash','%s/reformatblastoutput.sh'%sourcedir,outputpath, blastdbs, sourcedir,str(args.bidirectionalblast)], preexec_fn=default_sigpipe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr= p.communicate()
    if stdout.strip()=='no blast hits':
        noblasthits=True
        print('Warning: no blast alignments produced from genome comparison(s); pipeline terminated')
    if noblasthits==False:
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished reformatting alignments')
        print('finished reformatting alignments')
        runsubprocess(['bash','%s/getseqlengths.sh'%sourcedir,outputpath,blasttype,fastadir,sourcedir])
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished getting sequence lengths')
        print('finished getting sequence lengths')

    if args.keep==0 or args.keep==1:
        runsubprocess(['rm -rf %s/splitfastas'%outputpath],shell=True)        

if args.sequences.name=='<stdin>':
    blasttype='pairwise'
    fastadir1='%s/splitfastas1'%outputpath
    fastadir2='%s/splitfastas2'%outputpath
    fastafiles1='%s/fastafilepaths1.tsv'%outputpath
    blastdbs1='%s/blastdbfilepaths1.tsv'%outputpath
    fastafiles2='%s/fastafilepaths2.tsv'%outputpath
    blastdbs2='%s/blastdbfilepaths2.tsv'%outputpath
    if os.path.exists(fastadir1):
        sys.exit('Error: %s directory already exists, delete directory and try again'%fastadir)
    if os.path.exists(fastadir2):
        sys.exit('Error: %s directory already exists, delete directory and try again'%fastadir)
    splitfastas(args.sequences1,fastadir1,fastafiles1,blastdbs1)
    splitfastas(args.sequences2,fastadir2,fastafiles2,blastdbs2)
    #check there is no overlap between fastas provided in -s1 and -s2
    s1dir=os.listdir('%s/splitfastas1'%outputpath)
    s2dir=os.listdir('%s/splitfastas2'%outputpath)
    s1fastas=[f for f in s1dir if f.endswith('.fasta')]
    s2fastas=[f for f in s2dir if f.endswith('.fasta')]
    overlap=len(set(s1fastas).intersection(set(s2fastas)))
    if overlap>0:
        parser.error('there must be no overlap between fasta identifiers contained within the fasta files provided using the -s1 and -s2 flags')
    runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,fastafiles1, str(args.threads),sourcedir])
    runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,fastafiles2, str(args.threads),sourcedir])
    laterruntime=runtime()
    #print(laterruntime-startruntime, 'runtime; finished creating blast databases')
    print('finished creating blast databases')
    p=subprocess.Popen(['bash','%s/runblast.sh'%sourcedir,outputpath, fastadir1, blastdbs2, str(args.evalue), str(args.wordsize), str(args.task),str(args.threads)], preexec_fn=default_sigpipe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr= p.communicate()
    #try:
    #    print('{} {}'.format(stdout.decode(), 'stdout'))
    #except:
    #    pass
    #try:
    #    print('{} {}'.format(stderr.decode(), 'stderr'))
    #except:
    #    pass
    if p.returncode!=0:
        sys.exit('unexpected error when running runblast.sh')
    p=subprocess.Popen(['bash','runblast.sh',outputpath, fastadir2, blastdbs1, str(args.evalue), str(args.wordsize), str(args.task),str(args.threads),str(args.bidirectionalblast)], preexec_fn=default_sigpipe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr= p.communicate()
    #try:
    #    print('{} {}'.format(stdout.decode(), 'stdout'))
    #except:
    #    pass
    #try:
    #    print('{} {}'.format(stderr.decode(), 'stderr'))
    #except:
    #    pass
    if p.returncode!=0:
        sys.exit('unexpected error when running runblast.sh')
        
    laterruntime=runtime()
    #print(laterruntime-startruntime, 'runtime; finished running blast')
    print('finished running blast')
    
    blastdbs='%s/blastdbfilepaths_combined.tsv'%outputpath
    runsubprocess(['cat %s %s | sort -k1,1V > %s'%(blastdbs1,blastdbs2,blastdbs)],shell=True)

    p=subprocess.Popen(['bash','%s/reformatblastoutput.sh'%sourcedir,outputpath, blastdbs, sourcedir,str(args.bidirectionalblast)], preexec_fn=default_sigpipe, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr= p.communicate()
    if stdout.strip()=='no blast hits':
        noblasthits=True
        print('Warning: no blast alignments produced from genome comparison(s); pipeline terminated')
    if noblasthits==False:
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished reformatting alignments')
        print('finished reformatting alignments')
        runsubprocess(['bash','%s/getseqlengths.sh'%sourcedir,outputpath,blasttype,fastadir1,fastadir2,sourcedir])
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished getting sequence lengths')
        print('finished getting sequence lengths')

    if args.keep==0 or args.keep==1:
        runsubprocess(['rm -rf %s/splitfastas1'%outputpath],shell=True)
        runsubprocess(['rm -rf %s/splitfastas2'%outputpath],shell=True)


if args.blastonly!=True and noblasthits==False:
    runsubprocess(['Rscript','%s/granges.R'%sourcedir,outputpath, str(args.threads), str(args.breakpoint),str(args.alnlenstats),str(args.boot),str(args.alnrankmethod),str(args.lengthfilter),str(args.trimmedalignments),str(args.disjointalignments),str(args.bidirectionalblast)])
    laterruntime=runtime()
    #print(laterruntime-startruntime, 'runtime; finished trimming alignments')
    print('finished trimming alignments')
    #get included/excluded samples
    runsubprocess(['python','%s/getincludedexcluded.py'%sourcedir,outputpath,blastdbs])
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

          linecount=int(runsubprocess(['wc -l < %s/included.txt'%outputpath],shell=True,printstdout=False))
          if linecount < 3:
              print('Warning: there are only %i samples: a tree cannot be constructed'%linecount)
          else:
              treeargs=['Rscript','%s/tree.R'%sourcedir,outputpath,str(args.threads),str(args.boot),dendrogram,phylogeny]
              treeargs.extend(args.distscore)
              runsubprocess(treeargs)
              laterruntime=runtime()
              #print(laterruntime-startruntime, 'runtime; finished plotting tree(s) using distance score(s): %s'%args.distscore)
              print('finished plotting tree(s) using distance score(s): %s'%args.distscore)

if args.keep==0:
    runsubprocess(['rm -rf %s/blast'%outputpath],shell=True)

#OLD CODE
# parser = argparse.ArgumentParser(description='Run pipeline scripts')
# parser.add_argument('-s','--sequences', help='Sequences, for all-vs-all pairwise comparison', required=False)
# parser.add_argument('-s1','--sequences1', help='First set of sequence(s), for pairwise comparison against second set', required=False)
# parser.add_argument('-s2','--sequences2', help='Second set of sequence(s), for pairwise comparison against first set', required=False)
# parser.add_argument('-o','--out', help='Output directory (required)', required=True)
# parser.add_argument('-e','--evalue', help='BLAST e-value cutoff (default: 1e-8)', default=1e-8, type=float) #1e-8 is used in ggdc web pipeline - see Meier-Kolthoff 2014
# parser.add_argument('-w','--wordsize', help='BLAST word size (default: 38)', default=38, type=int) #38 is used in ggdc web pipleine?
# parser.add_argument('-t','--threads', help='Number of threads to use (default: 1)', default=1, type=int)
# parser.add_argument('-b','--boot', help='Number of bootstraps to run (default: no bootstrapping)', default=0, type=positiveint)
# parser.add_argument('-l','--lengthfilter', help='Length threshold (in basepairs) used to filter alignments prior to calculating breakpoint distance (default: 800)', default=800, type=positiveint)
# parser.add_argument('-d','--distscore', help='Distance score to use to construct tree (can specify multiple parameters; default: DistanceScore_d8 DistanceScore_d9)', nargs='+', default=['DistanceScore_d8', 'DistanceScore_d9'], choices=['DistanceScore_d0','DistanceScore_d4','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9'], metavar='',type=str)
# parser.add_argument('-m','--treemethod', help='Tree building method (can specify dendrogram and/or phylogeny, or none (i.e. no tree will be plotted); default: dendrogram)', nargs='+', default=['dendrogram'], choices=['dendrogram','phylogeny','none'], metavar='',type=str)
# parser.add_argument('-r','--alnrankmethod', help='Parameter used for selecting best alignment; default: bitscore)', default='bitscore', choices=['bitscore','alnlen','pid'], metavar='',type=str)
# parser.add_argument('--breakpoint', action='store_true', help='Calculate breakpoint statistics (default: do not calculate)')
# parser.add_argument('--alnlenstats', action='store_true', help='Calculate alignment length distribution statistics (default: do not calculate)')
# parser.add_argument('--trimmedalignments', action='store_true', help='Write to file trimmed alignments for each sample (default: do not output)')

#!/usr/bin/env python
import argparse, os, sys, time, gzip
from argparse import RawTextHelpFormatter

sourcedir=os.path.dirname(os.path.abspath(__file__))
cwdir=os.getcwd()
sys.path.append(sourcedir)

from pythonmods import runsubprocess,splitfastas

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

def fastainput(x):
    if os.path.exists(x): #file,directory,or (unbroken) symbolic link
        return str(x)
    else:
        return sys.stdin

parser = argparse.ArgumentParser(description="ATCG: Alignment Based Tool for Comparative Genomics",add_help=False,formatter_class=RawTextHelpFormatter)
#Help options
help_group = parser.add_argument_group('Help')
help_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
#Input options                                                               
input_group = parser.add_argument_group('Input')
input_group.add_argument('-s','--sequences', help='Sequences, for all-vs-all pairwise comparison (required if -1 and -2 flags not provided)', required=False,type=fastainput)
input_group.add_argument('-1','--sequences1', help='First set of sequence(s), for pairwise comparison against second set (required if -s flag not provided)', required=False,type=fastainput)
input_group.add_argument('-2','--sequences2', help='Second set of sequence(s), for pairwise comparison against first set (required if -s flag not provided)', required=False,type=fastainput)
#Output options                                               
output_group = parser.add_argument_group('Output')
output_group.add_argument('-o','--out', help='Output directory (required)', required=True)
output_group.add_argument('-d','--distscore', help='Distance score to use to construct tree; can specify multiple parameters (default: DistanceScore_d8 DistanceScore_d9)', nargs='+', default=['DistanceScore_d8', 'DistanceScore_d9'], choices=['DistanceScore_d0','DistanceScore_d4','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9'], metavar='', type=str)
output_group.add_argument('-m','--treemethod', help='Tree building method; can specify dendrogram and/or phylogeny, or none (i.e. no tree will be plotted); (default: none)', nargs='+', default=['none'], choices=['dendrogram','phylogeny','none'], metavar='', type=str)
output_group.add_argument('--breakpoint', action='store_true', help='If flag is provided, calculate breakpoint statistics (default: do not calculate)')
output_group.add_argument('--alnlenstats', action='store_true', help='If flag is provided, calculate alignment length distribution statistics (default: do not calculate)')
output_group.add_argument('--bestblastalignments', action='store_true', help='If flag is provided, write to file best blast alignments for each sample (default: do not output)')
output_group.add_argument('--nonoverlappingalignments', action='store_true', help='If flag is provided, write to file best alignments with non-overlapping query/subject ranges (only the intersection of query and subject alignment sets is given) (default: do not output)')
output_group.add_argument('--trimmedalignments', action='store_true', help='If flag is provided, write to file trimmed alignments for each sample (default: do not output)')
output_group.add_argument('--statsfromtrimmed', action='store_true', help='If flag is provided, statistics will be calculated from the trimmed alignments rather than from the non-verlapping alignments (default: calculate stats from non-overlapping alignments)') #DEPRECATED - recommended to use the default
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
output_group.add_argument('--keep', default=0, choices=[0,1,2], type=int, help='Level of file retention (default: 0)\n  '
                          ' 0 = do not retain blast database or blast table output,\n'
                          '   1 = do not retain blast database output,\n'
                          '   2 = retain all intermediate output')

args = parser.parse_args()
outputpath=os.path.relpath(args.out, cwdir)
startruntime=runtime()

#check sequence input flags used correctly
if args.sequences==None and args.sequences1==None and args.sequences2==None:
    parser.error('as input, you must either provide --sequences or both --sequences1 and --sequences2')
if args.sequences!=None and args.sequences1!=None and args.sequences2!=None:
    parser.error('as input, you must either provide --sequences or both --sequences1 and --sequences2')
if (args.sequences1==None and args.sequences2!=None) or (args.sequences1!=None and args.sequences2==None):
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

#check if fasta input is stdin, file, or directory; if file or stdin, open filehandle
if args.sequences!=None:
    #-s flag provided
    try:
        os.stat(args.sequences) #check file path exists
    except:
        sys.exit('Error: -s flag filepath not found or broken symlink')
    if os.path.isdir(args.sequences):
        fastafileinput='dir'
    elif os.path.isfile(args.sequences):
        fastafileinput='file'
        if args.sequences.endswith('.gz'):
            args.sequences=gzip.open(args.sequences,'rt')
        else:
            args.sequences=open(args.sequences,'r')
    else:
        fastafileinput='stdin'
        try:
            args.sequences=open(args.sequences,'r')
        except:
            sys.exit('Error: failed to open fasta file provided to -s flag; %s not a valid filepath?'%args.sequences)
else:
    #-1/-2 flags provided
    try:
        os.stat(args.sequences1)
    except:
        sys.exit('Error: -1 flag filepath not found or broken symlink')
    if os.path.isdir(args.sequences1):
        fastafile1input='dir'
    elif os.path.isfile(args.sequences1):
        fastafile1input='file'
        if args.sequences1.endswith('.gz'):
            args.sequences1=gzip.open(args.sequences1,'rt')
        else:
            args.sequences1=open(args.sequences1,'r')
    else:
        fastafile1input='stdin'
        try:
            args.sequences1=open(args.sequences1,'r')
        except:
            sys.exit('Error: failed to open fasta file provided to -1 flag; %s not a valid filepath?'%args.sequences1)
    #
    try:
        os.stat(args.sequences2)
    except:
        sys.exit('Error: -2 flag filepath not found or broken symlink')
    if os.path.isdir(args.sequences2):
        fastafile2input='dir'
    elif os.path.isfile(args.sequences2):
        fastafile2input='file'
        if args.sequences2.endswith('.gz'):
            args.sequences2=gzip.open(args.sequences2,'rt')
        else:
            args.sequences2=open(args.sequences2,'r')
    else:
        fastafile2input='stdin'
        try:
            args.sequences2=open(args.sequences2,'r')
        except:
            sys.exit('Error: failed to open fasta file provided to -2 flag; %s not a valid filepath?'%args.sequences2)
    if (fastafile1input=='dir' and fastafile2input!='dir') or (fastafile2input=='dir' and fastafile1input!='dir'):
        if fastafile1input!='dir':
            args.sequences1.close()
        if fastafile2input!='dir':
            args.sequences2.close()
        sys.exit('Error: fasta file input provided to flags -1 and -2 must be either both files or both directories (files can be provided as stdin)')



noblasthits=False

if os.path.exists(outputpath):
    sys.exit('Error: %s output directory already exists, delete directory and try again'%outputpath)

if args.sequences!=None:
    blasttype='allvallpairwise'
    blastdbdir='%s/blastdbs'%outputpath
    filepathinfo='%s/filepathinfo.tsv'%outputpath
    subjectsamples='%s/allsubjects.txt'%outputpath
    if fastafileinput=='file' or fastafileinput=='stdin':
        splitfastas(args.sequences,blastdbdir,filepathinfo,'%s/seqlengths.tsv'%outputpath)
        runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,filepathinfo,str(args.threads),sourcedir])
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished creating blast databases')
        print('finished creating blast databases')
        runsubprocess(['bash','%s/runblast.sh'%sourcedir,outputpath, blastdbdir, filepathinfo, str(args.evalue), str(args.wordsize), str(args.task),str(args.threads),str(args.bidirectionalblast),blasttype],preexec_fn='sigpipefix')
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished running blast')
        print('finished running blast')
    else:
        runsubprocess(['python','%s/getdirpaths.py'%sourcedir,args.sequences,filepathinfo,blastdbdir,blasttype])
        runsubprocess(['python','%s/getseqlengths.py'%sourcedir,'%s/seqlengths.tsv'%outputpath,filepathinfo])
        runsubprocess(['bash','%s/makeblastdbs_dirinput.sh'%sourcedir,filepathinfo,str(args.threads),sourcedir])
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished creating blast databases')
        print('finished creating blast databases')
        runsubprocess(['bash','%s/runblast_dirinput.sh'%sourcedir,outputpath,sourcedir, filepathinfo, filepathinfo,str(args.evalue), str(args.wordsize), str(args.task),str(args.threads),str(args.bidirectionalblast),blasttype],preexec_fn='sigpipefix')
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished running blast')
        print('finished running blast')

    noblasthits=runsubprocess(['bash','%s/reformatblastoutput.sh'%sourcedir,outputpath, subjectsamples, sourcedir,str(args.bidirectionalblast)],preexec_fn='sigpipefix',printstdout=False)
    if str(noblasthits)=='no blast hits':
        print('Warning: no blast alignments produced from genome comparison(s); pipeline terminated')
        noblasthits=True
    else:
        noblasthits=False
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished reformatting alignments')
        print('finished reformatting alignments')

    if args.keep==0 or args.keep==1:
        runsubprocess(['rm -rf %s/blastdbs'%outputpath],shell=True)

if args.sequences==None:
    blasttype='pairwise'
    blastdbdir1='%s/blastdbs1'%outputpath
    blastdbdir2='%s/blastdbs2'%outputpath
    filepathinfo1='%s/filepathinfo1.tsv'%outputpath
    filepathinfo2='%s/filepathinfo2.tsv'%outputpath
    filepathinfo='%s/filepathinfo.tsv'%outputpath
    subjectsamples='%s/allsubjects.txt'%outputpath
    if fastafile1input=='file' or fastafile1input=='stdin':
        splitfastas(args.sequences1,blastdbdir1,filepathinfo1,'%s/seqlengths1.tsv'%outputpath)
        splitfastas(args.sequences2,blastdbdir2,filepathinfo2,'%s/seqlengths2.tsv'%outputpath)
        #check there is no overlap between fastas provided in -s1 and -s2
        s1dir=os.listdir('%s/blastdbs1'%outputpath)
        s2dir=os.listdir('%s/blastdbs2'%outputpath)
        s1fastas=[f for f in s1dir if f.endswith('.fasta')]
        s2fastas=[f for f in s2dir if f.endswith('.fasta')]
        overlap=len(set(s1fastas).intersection(set(s2fastas)))
        if overlap>0:
            parser.error('there must be no overlap between fasta identifiers contained within the fasta files provided using the -s1 and -s2 flags')
        runsubprocess(['cat %s %s > %s'%('%s/seqlengths1.tsv'%outputpath,'%s/seqlengths2.tsv'%outputpath,'%s/seqlengths.tsv'%outputpath)],shell=True)
        runsubprocess(['cat %s %s > %s'%(filepathinfo1,filepathinfo2,filepathinfo)],shell=True)
        runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,filepathinfo, str(args.threads),sourcedir])
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished creating blast databases')
        print('finished creating blast databases')
        runsubprocess(['bash','%s/runblast.sh'%sourcedir,outputpath, blastdbdir1, filepathinfo2, str(args.evalue), str(args.wordsize), str(args.task),str(args.threads),str(args.bidirectionalblast),blasttype],preexec_fn='sigpipefix')
        if str(args.bidirectionalblast)=='True':
            runsubprocess(['bash','%s/runblast.sh'%sourcedir,outputpath, blastdbdir2, filepathinfo1, str(args.evalue), str(args.wordsize), str(args.task),str(args.threads),str(args.bidirectionalblast),'pairwiserun2'],preexec_fn='sigpipefix')
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished running blast')
        print('finished running blast')
    else:
        runsubprocess(['python','%s/getdirpaths.py'%sourcedir,args.sequences1,filepathinfo1,blastdbdir1,blasttype])
        runsubprocess(['python','%s/getdirpaths.py'%sourcedir,args.sequences2,filepathinfo2,blastdbdir2,blasttype])
        runsubprocess(['python','%s/getseqlengths.py'%sourcedir,'%s/seqlengths1.tsv'%outputpath,filepathinfo1])
        runsubprocess(['python','%s/getseqlengths.py'%sourcedir,'%s/seqlengths2.tsv'%outputpath,filepathinfo2])
        runsubprocess(['cat %s %s > %s'%('%s/seqlengths1.tsv'%outputpath,'%s/seqlengths2.tsv'%outputpath,'%s/seqlengths.tsv'%outputpath)],shell=True)
        runsubprocess(['cat %s %s > %s'%(filepathinfo1,filepathinfo2,filepathinfo)],shell=True)
        runsubprocess(['bash','%s/makeblastdbs_dirinput.sh'%sourcedir,filepathinfo, str(args.threads),sourcedir])
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished creating blast databases')
        print('finished creating blast databases')
        runsubprocess(['bash','%s/runblast_dirinput.sh'%sourcedir,outputpath,sourcedir,filepathinfo1, filepathinfo2,str(args.evalue), str(args.wordsize), str(args.task),str(args.threads),str(args.bidirectionalblast),blasttype],preexec_fn='sigpipefix')
        if str(args.bidirectionalblast)=='True':
            runsubprocess(['bash','%s/runblast_dirinput.sh'%sourcedir,outputpath,sourcedir,filepathinfo2, filepathinfo1,str(args.evalue), str(args.wordsize), str(args.task),str(args.threads),str(args.bidirectionalblast),'pairwiserun2'],preexec_fn='sigpipefix')
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished running blast')
        print('finished running blast')


    noblasthits=runsubprocess(['bash','%s/reformatblastoutput.sh'%sourcedir,outputpath, subjectsamples, sourcedir,str(args.bidirectionalblast)],preexec_fn='sigpipefix',printstdout=False)
    if str(noblasthits)=='no blast hits':
        print('Warning: no blast alignments produced from genome comparison(s); pipeline terminated')
        noblasthits=True
    else:
        noblasthits=False
        laterruntime=runtime()
        #print(laterruntime-startruntime, 'runtime; finished reformatting alignments')
        print('finished reformatting alignments')

    if args.keep==0 or args.keep==1:
        runsubprocess(['rm -rf %s/blastdbs1'%outputpath],shell=True)
        runsubprocess(['rm -rf %s/blastdbs2'%outputpath],shell=True)


if args.blastonly!=True and noblasthits==False:
    runsubprocess(['Rscript','%s/granges.R'%sourcedir,outputpath, str(args.threads), str(args.breakpoint),str(args.alnlenstats),str(args.boot),str(args.alnrankmethod),str(args.lengthfilter),str(args.bestblastalignments),str(args.nonoverlappingalignments),str(args.trimmedalignments),str(args.bidirectionalblast),str(args.statsfromtrimmed)])
    laterruntime=runtime()
    #print(laterruntime-startruntime, 'runtime; finished parsing alignments')
    print('finished parsing alignments')
    #get included/excluded samples
    runsubprocess(['python','%s/getincludedexcluded.py'%sourcedir,outputpath,filepathinfo])
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

if args.keep==0 or args.keep==1:
    runsubprocess(['rm %s/includedsubjects.txt'%outputpath],shell=True)
    runsubprocess(['rm %s/allsubjects.txt'%outputpath],shell=True)
    if args.sequences==None:
        runsubprocess(['rm %s/filepathinfo.tsv'%outputpath],shell=True) #filepathinfo1.tsv/filepathinfo2.tsv retained though
        runsubprocess(['rm %s/seqlengths1.tsv'%outputpath],shell=True)  #seqlengths.tsv retained
        runsubprocess(['rm %s/seqlengths2.tsv'%outputpath],shell=True)


print('Finished running ATCG!')

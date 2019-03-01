#!/usr/bin/env python                                                                                                                                                                                      
import argparse, os, sys, subprocess
sourcedir=os.path.dirname(os.path.abspath(__file__))
cwdir=os.getcwd()
sys.path.append(sourcedir)

from pythonmods import runsubprocess

###gffparsing.py takes either multigenbank file or concatenated genbank files (as stdout from concatenategenbank.sh); outputs per-genome gffs using gffparsing2.py

###fixprokka.sh takes either a multigff file from prokka or concatenated prokka gff files (as stdout from concatenateprokka.sh); outputs per-genome gffs

parser = argparse.ArgumentParser(description="ATCG: Alignment Based Tool for Comparative Genomics; get feature annotation files in correct format for visualisation.py",add_help=False)
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
parser.add_argument('-t', '--annotationtype', help='The type of annotation file that requires conversion to correct format (required)',choices=['prokka','genbank'],type=str,required=True)
parser.add_argument('-i', '--inputpath', help='The input directory (containing annotation files) or annotation file to be converted to correct format (required)',required=True)
parser.add_argument('-o', '--outdir', help='The output directory (required)',required=True)
parser.add_argument('-s', '--seqnames', help='A file containing the sequence names associated with the annotation file(s) in the first column (required if annotationtype is prokka)',required=False)

args = parser.parse_args()
outputpath=os.path.relpath(args.outdir, cwdir)

if args.seqnames==None:
    if args.annotationtype=='prokka':
        print('Error: if using prokka annotation file(s) as input, a file containing the associated sequence names (the original names, with no changes introduced by prokka) must be provided')
        sys.exit()

runsubprocess(['mkdir -p %s'%outputpath],shell=True)

#check if input is file or directory 
if os.path.isfile(args.inputpath):
    inputpathtype='file'                                                                                 
elif os.path.isdir(args.inputpath):
    inputpathtype='directory'
else:
    print('Error: %s is not a file or directory'%args.inputpath)                                                                                              
    sys.exit()

if args.annotationtype=='prokka':
    if inputpathtype=='directory':
        runsubprocess(['bash %s/concatenateprokka.sh %s | python %s/fixprokkagff.py %s %s %s'%(sourcedir,str(args.inputpath),sourcedir,str(args.seqnames),outputpath,inputpathtype)],shell=True)
    else:
        runsubprocess(['python', '%s/fixprokkagff.py'%sourcedir, str(args.seqnames),outputpath,inputpathtype,str(args.inputpath)])
else:
    if inputpathtype=='directory':
        runsubprocess(['bash %s/concatenategbk.sh %s | python %s/gffparsing.py %s | python %s/gffparsing2.py %s'%(sourcedir,str(args.inputpath),sourcedir,inputpathtype,sourcedir,outputpath)],shell=True)
    else:
        runsubprocess(['python %s/gffparsing.py %s | python %s/gffparsing2.py %s'%(sourcedir,inputpathtype,sourcedir,outputpath)],shell=True)   

        

#OLD NOTES
#genbanktogff; need to incorporate gffparsing as function - need to handle multigenbank (not necessary)

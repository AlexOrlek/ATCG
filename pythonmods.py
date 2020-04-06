#Python modules

###wrapper for subprocess command
def runsubprocess(args,verbose=False,shell=False,polling=False,printstdout=True,preexec_fn=None):
    """takes a subprocess argument list and runs Popen/communicate or Popen/poll() (if polling=True); if verbose=True, processname (string giving command call) is printed to screen (processname is always printed if a process results in error); errors are handled at multiple levels i.e. subthread error handling; fuction can be used fruitfully (returns stdout)"""
    #function setup
    import subprocess,sys,signal
    try:
        import thread
    except:
        import _thread
    def subprocess_setup(): #see: https://github.com/vsbuffalo/devnotes/wiki/Python-and-SIGPIPE
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    if preexec_fn=='sigpipefix':
        preexec_fn=subprocess_setup
    if shell==True:
        processname=args[0]
        processname=processname[0].split()
        processname=(" ".join(a for a in processname))
    else:
        processname=(" ".join(a for a in args))
    #
    if verbose==True:
        print('{} {}'.format(processname, 'processname'))
    try:
        if polling==True:
            p=subprocess.Popen(args, stdout=subprocess.PIPE,shell=shell,preexec_fn=preexec_fn)
            while True:
                stdout=p.stdout.readline()
                if p.poll() is not None:
                    break
                if stdout: #if stdout not empty...
                    if printstdout==True:
                        print('{}'.format(stdout.decode().strip()))
        else:
            p=subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=shell,preexec_fn=preexec_fn)
            stdout, stderr= p.communicate()
            if stdout:
                if printstdout==True:
                    print('{}'.format(stdout.decode()))
            if stderr:
                try: #want to output to stderr stream
                    if (sys.version_info > (3, 0)):
                        print('{}'.format(stderr.decode()),file=sys.stderr) #Python3
                    else:
                        print>>sys.stderr,stderr  #Python2
                except: #if above code block fails for some reason, print stderr (to stdout)
                    print('{}'.format(stderr.decode()))

        if p.returncode==0:
            if verbose==True:
                print('{} {}'.format(processname, 'code has run successfully'))
        else:
            sys.exit() #triggers except below
    except:
        print('{} {}'.format(processname, '#this pipeline step produced error'))
        print('unexpected error; exiting')
        sys.exit()
        
    if p.returncode!=0:
        print('unexpected error; exiting')
        try:
            thread.interrupt_main()
        except:
            _thread.interrupt_main()
    else:
        if stdout:
            return stdout.decode().strip()



###wrapper for splitting input fasta by sample
def splitfastas(infile,fastadir,fastafilepaths,blastdbfilepaths):
    """takes input fastafile from filepath or sys.stdin; splits by sample and writes to outdir; also writes fastafilepaths and blastdbfilepaths"""
    from Bio import SeqIO
    import re,os
    from pythonmods import runsubprocess
    #first create splitfasta output directory
    runsubprocess(['mkdir -p %s'%fastadir],shell=True)
    #parse fasta and store in recorddict
    recorddict={}
    for recordid,recordseq in SeqIO.FastaIO.SimpleFastaParser(infile):
        #remove description from fasta id if present
        newfastaheader=re.sub(r'(\S+)(?: .*)?',r'\1',recordid)
        newfastaheader=newfastaheader.strip()
        recordid=newfastaheader
        #get sample name (fasta header format should be sample or sample|contig)
        sample=re.match(r'^([^\|]*).*',newfastaheader)
        sample=sample.group(1)
        #write to dict
        if sample not in recorddict:
            recorddict[sample]=[]
        recorddict[sample].append((recordid,recordseq))
    #write records to splitfastas directory, split by sample
    f2=open(fastafilepaths,'w')
    f3=open(blastdbfilepaths,'w')
    for sample in recorddict.keys():
        fastafilepath='%s/%s.fasta'%(fastadir,sample)
        blastdbfilepath=os.path.splitext(fastafilepath)[0]
        blastdbfilepath='%s_db'%blastdbfilepath
        f2.write('%s\t%s\n'%(sample,fastafilepath))
        f3.write('%s\t%s\n'%(sample,blastdbfilepath))
        with open(fastafilepath,'w') as output_handle:
            for recordid,recordseq in recorddict[sample]:
                output_handle.write(">%s\n%s\n" % (recordid, recordseq))
    f2.close()
    f3.close()


###wrappers for blast commands
def makeBLASTdb(fastafile, databasename, dbtype, parse_seqids=False): #dbtype can be 'nucl' or 'prot'
    """takes fastafile filepath, databasename filepath and dbtype args"""
    import subprocess
    if parse_seqids==False:
        cmdArgs=['makeblastdb', '-dbtype', dbtype, '-in', fastafile, '-out', databasename]
    else:
        cmdArgs=['makeblastdb', '-dbtype', dbtype, '-in', fastafile, '-out', databasename, '-parse_seqids']
    subprocess.call(cmdArgs)
            


def runblastn(query, database, blastoutput, evalue=str(10), outfmt='custom qcov', task='blastn', num_threads=str(1), max_target_seqs=str(500), max_hsps=False, perc_identity=False, qcov_hsp_perc=False, culling_limit=False, word_size=False): #default evalue is 10; default wordsize is 11 for blastn - just use -tasks parameter which also changes gap settings; default task with no -task parameter set is megablast; N.b use cutstom qcov as standard for blast-based typing and plasmidpipeline, otherwise use outfmt 6 or specify custom outfmt
    import subprocess
    evalue=str(evalue)
    if outfmt=='custom qcov':
        outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen'
    cmdArgs=['blastn', '-query', query, '-db', database, '-out', blastoutput, '-evalue', evalue, '-outfmt', outfmt, '-task', task, '-num_threads', num_threads, '-max_target_seqs', max_target_seqs]
    if max_hsps!=False:
        cmdArgs.extend(['-max_hsps', max_hsps])
    if perc_identity!=False:
        cmdArgs.extend(['-perc_identity', perc_identity])
    if qcov_hsp_perc!=False:
        cmdArgs.extend(['-qcov_hsp_perc', qcov_hsp_perc])
    if culling_limit!=False:
        cmdArgs.extend(['-culling_limit', culling_limit])
    if word_size!=False:
        cmdArgs.extend(['-word_size', word_size])
    subprocess.call(cmdArgs)




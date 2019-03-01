import sys,os
data = sys.stdin.read().strip().split('\n')

for line in data:
    line=line.split('\t')
    seqname=line[0].split()[0]
    seqname=seqname.split('|')
    if len(seqname)>2:
        seqname=seqname[0:2]
    seqname='|'.join(seqname)
    seqlen=line[1]
    line='\t'.join([seqname,seqlen])
    print(line)
    

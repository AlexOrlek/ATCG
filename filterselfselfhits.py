import sys

infile=sys.stdin.read()

infile=infile.strip().split('\n')
for line in infile:
    data=line.split('\t')
    query=data[0].split('|')[0]
    subject=data[1].split('|')[0]
    if query != subject:
        print('%s'%'\t'.join(data))

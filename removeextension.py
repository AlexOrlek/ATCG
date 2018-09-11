import sys,os
data = sys.stdin.read().split()

for line in data:
    line=line.strip()
    line=os.path.splitext(line)[0]
    print(line)

#!/usr/bin/env python3
'''
wrapper for last program
run last  with BlastTab+ output and return mgblast like formated output
last BlastTab+ output column order:
1 query name
2 reference name
3 percent identity
4 alignment length
5 mismatches
6 gap opens
7 query start
8 query end
9 reference start
10 reference end
11 e-value
12 bitscore
13 length of query
14 length of reference sequence
(accordin lastal manual - more column may be added in future)

Needed mgblast order:

qseqid   1 -> 1
qlen     2 -> 13
qstart   3 -> 7
qend     4 -> 8
sseqid   5 -> 2
slen     6 -> 14
sstart   7 -> 9
send     8 -> 10
pident   9 -> 3
bitscore 10-> 12
evalue   11-> 11
sstrand  12-> must be evaluated!

'''
import subprocess
import sys
last_command = " ".join(["lastal"] + sys.argv[1:])
p = subprocess.Popen(last_command, shell=True, stdout=subprocess.PIPE)
for j in p.stdout:
    line = j.decode()
    if line[0] != "#":
        items = line.split("\t")
        strand = "+" if int(items[6]) < int(items[7]) else "-"
        out = "\t".join([items[i - 1]
                         for i in [1, 13, 7, 8, 2, 14, 9, 10, 3, 12, 11]
                         ]) + "\t" + strand + "\n"
        print(out, end="")

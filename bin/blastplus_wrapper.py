#!/usr/bin/env python3

import subprocess
import sys
blast_command = " ".join(sys.argv[1:])
print(blast_command, file=sys.stderr)
p = subprocess.Popen(blast_command, shell=True, stdout=subprocess.PIPE)
for j in p.stdout:
    line = j.decode()
    if line[0] != "#":
        items = line.strip().split("\t")
        # skip self hits:
        if items[0] == items[4]:
            continue
        items[11] = "+" if items[11] == "plus" else "-"
        out = "\t".join(items)
        print(out)

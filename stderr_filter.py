#!/usr/bin/env python3
'''
Purpose of this script is to filters some massages output
stderr to prevent galaxy to raise the error
'''
import sys

string_to_detect = [
    'Karlin-Altschul parameters',
    'slippage may introduce errors',
    'Examining 5 or more matches is recommended',
    'DeprecationWarning: The binary mode of fromstring is deprecated',
]

string_to_remove = [
    ('error', 'errour'),
    ('warning', 'alert')
]
input_file = sys.argv[1]

with open(input_file) as f:
    for line in f:
        for s in string_to_detect:
            if s in line:
                new_line = "--" + line.lower()
                for r in string_to_remove:
                    new_line = new_line.replace(r[0], r[1])
                line = new_line
        print("parsed line:", line, file=sys.stderr)

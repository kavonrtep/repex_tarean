#!/usr/bin/env python3
import os
import hashlib

from itertools import chain



def md5checksum(filename, fail_if_missing=True):
    try:
        md5 = hashlib.md5()
        with open(filename, "rb") as f:
            for i in iter(lambda: f.read(4096), b""):
                md5.update(i)
    except FileNotFoundError as e:
        if not fail_if_missing:
            return "Not calculated!!!!  File {} is missing".format(filename)
        else:
            raise e

    return md5.hexdigest()


class FilePath(str):
    '''
    Extension of str - it just contain additional atribute showing that the string is alsp path to file
    '''

    def __new__(cls, string):
        obj = super(FilePath, cls).__new__(cls, string)
        obj.filepath = True
        return obj

    def relative(self, start):
        ''' return path relative to start'''
        return os.path.relpath(self, start)


def save_as_table(d, path, header=None, relative=True):
    ''' takes list of dictionaries and save csv file
    define header if you want to use specific order!
    '''
    pathdir = os.path.dirname(path)
    if not header:
        
        all_keys = [i.keys() for i in d]
        header = set(chain(*all_keys))
        print("header: ---------", header)
    with open(path, 'w') as f:
        f.write("\t".join(header))
        f.write("\n")
        for i in d:
            istr = []
            for key in header:
                if isinstance(i[key], FilePath):
                    if relative:
                        istr.append('"' + str(i[key].relative(pathdir)) + '"')
                    else:
                        istr.append('"' + str(i[key]) + '"')
                else:
                    if isinstance(i[key], str):
                        istr.append('"' + str(i[key] + '"'))
                    else:
                        istr.append(str(i[key]))

            f.write("\t".join(istr))
            f.write("\n")


def export_tandem_consensus(clusters_info, path, rank=1, n=1):
    ''' export tr consensu to file'''
    print("exporting fasta files")
    print(clusters_info)
    s = None
    with open(path, 'w') as f:
        for cl in clusters_info:
            print(cl)
            print(dir(cl))
            if cl.TR_consensus and rank == cl.tandem_rank:
                s = ">CL{index}_TR_{n}_x_{L}nt\n{sequence}\n".format(
                    index=cl.index,
                    n=n,
                    L=cl.TR_monomer_length,
                    sequence=n * cl.TR_consensus.replace('<pre>', ''))
                f.write(s)
    if s:
        return path
    else:
        return None


def file_len(filename):
    '''count number of lines in file'''
    with open(filename) as f:
        i = 0
        for i in f:
            i += i
    return i

def go2line(f, L):
    ''' find line L in file object f '''
    f.seek(0)
    if L == 0:
        return
    i = 0
    pos = f.tell()
    for line in f:
        i += 1
        if i == L:
            f.seek(pos)
            return
        else:
            pos = pos + len(line)

def format_query(x):
    '''
    make list for query in format ("x","y","x",...)
    '''
    out = '("'+ '","'.join(
        map(str, x)
    ) + '")'
    return out

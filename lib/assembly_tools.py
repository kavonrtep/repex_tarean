#!/usr/bin/env python3
import sys
import logging
import subprocess
import os
import tempfile
import config
import shutil
import itertools
import pickle
import shlex
from lib.parallel.parallel import parallel2 as parallel
from lib.parallel.parallel import get_max_proc

REQUIRED_VERSION = (3, 4)
MAX_BUFFER_SIZE = 100000
if sys.version_info < REQUIRED_VERSION:
    raise Exception("\n\npython 3.4 or higher is required!\n")
LOGGER = logging.getLogger(__name__)
MAX_FASTA_IN_DIRECTORY = 1000

def assembly(sequences, hitsort, clusters_info, assembly_dir,
             contigs_file, min_size_of_cluster_for_assembly = 0):
    '''
    Runs assembly on sequences (SequenceSet). Assembly is
    performed on each cluster separatelly, clusters are taken
    from hitsort(Graph)
    Cluster_listing - list of Clusters
    if cluster.tandem_rank is 1 or 2 - no assembly is performed!!
    '''

    # iterate over large clusters, assembly is performed on sequences stored in
    # cluster_little[index].fasta_file_full - for annotated clusters
    # sequences of small clusters are retrieved from database
    fasta_seqs = [
        i.fasta_file_full
        for i in clusters_info
        if not i.tandem_rank in config.SKIP_CAP3_ASSEMBLY_TANDEM_RANKS
    ]
    prefixes = ["CL{}".format(i.index)
                for i in clusters_info
                if not i.tandem_rank in config.SKIP_CAP3_ASSEMBLY_TANDEM_RANKS]
    LOGGER.info("Number of clusters for assembly: {}".format(hitsort.number_of_clusters))
    LOGGER.info("Assembling large clusters")
    assembly_files = parallel(cap3worker, fasta_seqs, prefixes)
    LOGGER.info("Large clusters assembled")
    # some clusters - tanem rank 1 assembled
    j = 0
    for cl in clusters_info:
        if cl.tandem_rank in config.SKIP_CAP3_ASSEMBLY_TANDEM_RANKS:
            cl.assembly_files = {i: None for i in config.CAP3_FILES_MAPPING}
            consensus_file = cl.dir + "/tarean_contigs.fasta"
            cl.assembly_files["{}.{}.contigs"] = consensus_file
            with open(cl.dir_tarean + "/tarean_contigs.fasta",
                      'r') as fin, open(consensus_file, 'w') as fout:
                for line in fin:
                    if line[0] == ">":
                        line = ">CL{}Contig{}".format(cl.index, line[1:])
                    fout.write(line)
        else:
            cl.assembly_files = assembly_files[j]
            j += 1
    # assembly of small clusters:
    # connection to files were results will be concatenated
    LOGGER.info("Assembly of small cluster - making tmp files")

    tmp_dir_root = tempfile.mkdtemp()
    tmp_subdir = tempfile.mkdtemp(dir=tmp_dir_root)
    nproc = get_max_proc()
    prefixes = [[] for i in range(nproc)]
    tmp_seq_files = [[] for i in range(nproc)]

    LOGGER.info("Assembly of small clusters - saving small cluster to tmp files")
    seq_dictionary = sequences.toDict()
    fasta_count = 0
    chunk_counter = itertools.cycle(range(nproc))
    for index in range(len(clusters_info) + 1, hitsort.number_of_clusters):
        chunk = next(chunk_counter)
        ids = hitsort.get_cluster_reads(index)
        if len(ids) < min_size_of_cluster_for_assembly:
            break
        prefixes[chunk].append("CL{}".format(index))
        fasta_count += 1
        if fasta_count > MAX_FASTA_IN_DIRECTORY:
            # create new subdir to keep number of files in directory low
            fasta_count = 1
            tmp_subdir = tempfile.mkdtemp(dir=tmp_dir_root)
        fasta_file_name = "{}/{}".format(tmp_subdir, index)
        write_seqDict_to_fasta(file_name=fasta_file_name, sequences=seq_dictionary, subset=ids)
        tmp_seq_files[chunk].append(fasta_file_name)
    del seq_dictionary
    LOGGER.info("Assembly of small clusters running")
    pickled_fparts_small_contigs = parallel(cap3worker_multiple, tmp_seq_files, prefixes)
    LOGGER.info("Assembly of small clusters done")
    # append to all results
    LOGGER.info("Assembly of small clusters - appending results")
    all_con = {}
    for fkey in config.CAP3_FILES_MAPPING:
        if config.CAP3_FILES_MAPPING[fkey]:
            all_con[fkey] = open("{}/{}".format(
                assembly_dir, config.CAP3_FILES_MAPPING[fkey]), 'w')
        else:
            all_con[fkey] = None
    for pickled_fparts_multiple in pickled_fparts_small_contigs:
        with open(pickled_fparts_multiple, 'rb') as fp:
            fparts_multiple = pickle.load(fp)
        os.unlink(pickled_fparts_multiple)
        for fparts in fparts_multiple:
            for fkey in fparts:
                if all_con[fkey]:
                    with open(fparts[fkey]) as f:
                        all_con[fkey].write(f.read())
                    os.unlink(fparts[fkey])
                else:
                    os.unlink(fparts[fkey])
    
    for fkey in all_con:
        if all_con[fkey]:
            all_con[fkey].close()
    LOGGER.info("Assembling - copying files to final location")
    # copy all contigs to one location: - contigs_file
    with open(contigs_file, 'w') as fout:
        for cl in clusters_info:
            # append contig from large clusters
            with open(cl.assembly_files["{}.{}.contigs"], 'r') as fin:
                for i in fin:
                    fout.write(i)
        # append contigs from small clusters
        small_aln_file = "{}/{}".format(
            assembly_dir, config.CAP3_FILES_MAPPING['{}.{}.aln'])
        # calculate profile on aln file of small clusters
        file_base_name = small_aln_file[:-4]
        os.system("align_parsing.pl -i {fn} -o {out}.info.fasta -p {out}.profile 2>&1".format(fn=small_aln_file, out=file_base_name))
        os.system("select_and_sort_contigs.pl {fn}.info.fasta 5 2>&1".format(fn=file_base_name))
        small_contig_file = file_base_name + ".info.fasta"
        with open(small_contig_file, 'r') as fin:
            for i in fin:
                fout.write(i)
    shutil.rmtree(tmp_dir_root)

def write_seqDict_to_fasta(file_name, sequences, subset):
    with open(file_name, 'w') as f:
        for i in subset:
            f.write(">{}\n{}\n".format(i, sequences[i]))




def cap3worker(seqfile, prefix="cap", cap3args=" -p 80 -o 40 "):
    prefix2 = "cap"
    cmd = "cap3 " + seqfile + cap3args + " -x " + prefix2
    with open(seqfile + "." + prefix2 + ".aln", "w") as aln:
        subprocess.check_call(shlex.split(cmd), shell=False, stdout=aln)
    # this generate three files
    files_dict = {}
    for fkey in config.CAP3_FILENAMES:
        fn = fkey.format(seqfile, prefix2)
        fn_tmp = "{}.tmp".format(fn)
        files_dict[fkey] = fn
        if config.CAP3_PATTERNS_REPLACE[fkey]:
            pattern, replacement = config.CAP3_PATTERNS_REPLACE[fkey]
            with open(fn, "r") as con_in, open(fn_tmp, 'w') as con_out:
                for line in con_in:
                    con_out.write(line.replace(pattern, replacement.format(
                        prefix)))
            os.rename(fn_tmp, fn)
    
    # make new meaningful names here

    for fkey in config.CAP3_FILES_GOODNAMES:
        config.CAP3_FILES_GOODNAMES[fkey]
        fn_goodname = os.path.dirname(files_dict[fkey]) + "/" + config.CAP3_FILES_GOODNAMES[fkey]
        os.rename(files_dict[fkey], fn_goodname)
        files_dict[fkey] = fn_goodname

    aln_file = files_dict["{}.{}.aln"]
    file_base_name = aln_file[:-4]
    os.system("align_parsing.pl -i {fn} -o {out}.info.fasta -p {out}.profile 2>&1".format(fn=aln_file, out=file_base_name))
    os.system("select_and_sort_contigs.pl {fn}.info.fasta 5".format(fn=file_base_name))
    # TODO -add new file to files_dict
    # replace simple fasta with info.fasta
    files_dict["{}.{}.contigs"] = file_base_name + ".info.fasta"
    return files_dict

def cap3worker_multiple(many_seqfile, many_prefixes, cap3args=" -p 80 -o 40 "):
    '''
    purpose of this script is to run multiple assemblies within single process
    avoiding running high number of short parallel subprocesses
    As starting subprocess for each cap3 assembly was very ineffective,
    all ap3 commands are written to single file and run using single subprocess
    command -
    '''
    cmd_file = tempfile.NamedTemporaryFile(mode="w",delete=False).name
    with open(cmd_file,'w') as cmdf:
        for seqfile, prefix in zip(many_seqfile, many_prefixes):
            cmd = "cap3 " + seqfile + cap3args + " -x " + prefix + " > " + seqfile + "." + prefix + ".aln\n"
            cmdf.write(cmd)
    os.system("sh "+ cmd_file)
    # collect results:
    files_dict_many = []
    for seqfile, prefix in zip(many_seqfile, many_prefixes):
        files_dict = {}
        for fkey in config.CAP3_FILENAMES:
            fn = fkey.format(seqfile, prefix)
            fn_tmp = "{}.tmp".format(fn)
            files_dict[fkey] = fn
            if config.CAP3_PATTERNS_REPLACE[fkey]:
                pattern, replacement = config.CAP3_PATTERNS_REPLACE[fkey]
                with open(fn, "r") as con_in, open(fn_tmp, 'w') as con_out:
                    for line in con_in:
                        con_out.write(line.replace(pattern, replacement.format(
                            prefix)))
                os.rename(fn_tmp, fn)
        files_dict_many.append(files_dict)
    # this is too large to be return directly - use picking
    f = tempfile.NamedTemporaryFile(delete=False)
    os.unlink(cmd_file)
    pickle.dump(files_dict_many,file=f)
    return f.name
    


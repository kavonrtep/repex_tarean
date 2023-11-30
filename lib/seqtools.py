#!/usr/bin/env python3
import logging
logger = logging.getLogger(__name__)
import itertools
import os
import sys
import random
import sqlite3
import subprocess
import shlex  # for command line arguments split
from collections import namedtuple, OrderedDict
from lib.utils import format_query
from lib.parallel.parallel import parallel2 as parallel
from lib.parallel.parallel import get_max_proc
REQUIRED_VERSION = (3, 4)
MAX_PRINT = 10

if sys.version_info < REQUIRED_VERSION:
    raise Exception("\n\npython 3.4 or higher is required!\n")

# additional functions



def _hitsort_worker(query, blastdb, output, options):

    # output from blast is parsed from stdout
    cmd = options.all2all_search_params.format(query = query, blastdb = blastdb)
    print(cmd)
    min_lcov, min_pid, min_ovl, min_scov, evalue_max = options.filtering_threshold
    pair2insert = ''
    signs ={'plus':'+', 'minus':'-'}
    with open(output, 'w') as f:
        with subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE) as p:
            for i in p.stdout:
                items = i.decode('utf-8').split()
                if options.filter_self_hits:
                    if items[4] == items[0]:
                        continue
                evalue = float(items[10])
                ovl_q = abs(float(items[2]) - float(items[3])) + 1
                ovl_h = abs(float(items[6]) - float(items[7])) + 1
                if (ovl_q >= min_ovl or ovl_h >= min_ovl) and float(items[8]) >= min_pid:
                    if float(items[1]) >= float(items[5]):
                        lcov = ovl_q * 100.0 / float(items[1])
                        scov = ovl_h * 100.0 / float(items[5])
                    else:
                        lcov = ovl_q * 100.0 / float(items[5])
                        scov = ovl_h * 100.0 / float(items[1])
                    # positive line:
                    if lcov >= min_lcov and scov >= min_scov and evalue_max > evalue:
                        pr = [items[0], items[4]]
                        # previous HSP
                        if pair2insert != "{0}\t{1}".format(pr[0], pr[1]):
                            pair2insert = "{0}\t{1}".format(pr[0], pr[1])
                            if items[11] in ['plus', 'minus']:
                                items[11] = signs[items[11]]
                            f.write("{0}\t{1}\t{2}\n".format(pair2insert, items[9], "\t".join(
                                items[k] for k in [1, 2, 3, 5, 6, 7, 8, 10, 11])))

def blast_with_filter(fasta_file_filter, blastdb):
    '''
    Return list of sequences query id which has similarity to filter
    It uses mgblast for search
    '''
    params = '-p 85 -W18 -UT -X40 -KT -JF  -F "m D" -v100000000 -b100000000 -D4 -C 30 -H 30'
    positive_hits = set()
    min_pid = 90
    min_ovl_perc = 90
    cmd = " ".join(['mgblast',
                    '-i', fasta_file_filter,
                    '-d', blastdb,
                    params
                ])
    with subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE) as p:
        for i in p.stdout:
            items = i.decode('utf-8').split()
            ovl_perc = (abs(float(items[6]) - float(items[7])) + 1) / float(items[5]) * 100
            pid = float(items[8])
            if pid > min_pid and ovl_perc > min_ovl_perc:
                positive_hits.add(items[4])
    return(positive_hits)



# TODO test task
# predefine blast params
def blast_worker(query, blastdb, output, params):
    if params['program'] in ['blastx', 'blastn']:
        default_columns = ' -outfmt "6 {}"'.format(params['output_columns'])
        cmd = "{} -query {} -db {} {} {}".format(
            params['program'], query, blastdb, params['args'], default_columns)
    elif params['program']=='diamond blastx':
        # diomand have slightly different format than blastx
        default_columns = ' --outfmt 6 {}'.format(params['output_columns'])
        cmd = "{} -q {} -d {} {} {}".format(
            params['program'], query, blastdb, params['args'], default_columns)
    print(cmd)
    preceeding_pair = ['', '']
    BlastValues = namedtuple("BlastValues", params['output_columns'])
    with open(output, 'wb') as f:
        with subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE) as p:
            for i in p.stdout:
                items = i.decode('utf-8').split()
                blast_values = BlastValues(*[f(i) for i, f in zip(items, params['column_types'])])
                #qseqid, sseqid, qlen, slen, length, ppos, bitscore = [
                #    f(i) for i, f in zip(items, params['column_types'])]
                if params['filter_function'](blast_values):
                    if preceeding_pair != [blast_values.qseqid, blast_values.sseqid]:
                        f.write(i)
                        preceeding_pair = [blast_values.qseqid, blast_values.sseqid]


class Sequence:

    def __init__(self, seq, name, paired=False):
        # the mode os seq storage can be changed later to make it more
        # memory efficient
        self._seq = bytes(str(seq), "ascii")
        self.name = str(name)

    @property
    def seq(self):
        return self._seq.decode("utf-8")

    @seq.setter
    def seq(self, value):
        self._seq = bytes(str(value), "ascii")

    def __str__(self):
        return "{0} : {1}".format(self.name, self.seq)

    @staticmethod
    def read_fasta(fasta_file_name):
        '''
        generator - reads sequences from fasta file
        return sequence one by one
        '''
        with open(fasta_file_name, 'r') as f:
            header = None
            seqstr = None
            for rawline in f:
                line = rawline.strip()
                if line == "":
                    continue
                if line[0] == ">":
                    if header and seqstr:
                        yield Sequence(seqstr, header)
                        # reset
                        seqstr = None
                        header = line[1:]
                    elif seqstr:
                        Warning("sequence was not preceeded by header")
                    else:
                        header = line[1:]
                else:
                    seqstr = line if not seqstr else seqstr + line
        # skip empty lines:
        if header and seqstr:
            yield Sequence(seqstr, header)
        return

    def write2fasta(self, file_object):
        file_object.write(">{0}\n{1}\n".format(self.name, self.seq))


class SequenceSet:

    def __init__(self, filename=None, source=None, sample_size=0, new=False, paired=False, prefix_length=0, rename=False, seed=123, fasta=None):
        '''
        filename: path to sqlite database, if none database is stored in memory
        source: path to fasta file, if none empty database is created
        sample_size: use only sample of sample_size from fasta file, if 0 all is used
        new: should be old database created?
        paired - paired reads
        prefix_length -int
        rename -- use int as names

        creates SequenceSet, either empty or from fasta file,
        it can be stored as dictionary or in shelve in file
        only sample of all sequences can be loaded. if sample_size = 0
        all sequences are loaded
        if source is give, new database is created - if filename exist - it will be deleted!
        '''
        #======================================================================
        # TODO name checking
        # TODO detect unique prefixes
        #======================================================================
        # type of storage
        random.seed(seed)

        self.source = source
        self.sample_size = sample_size
        self.paired = paired
        self.filename = filename
        self.prefix_length = prefix_length
        self.prefix_codes = {}
        self.rename = rename
        self.fasta = fasta
        self._length = None
        self._original_length = None
        self.blastdb = None
        self.blastdb_legacy = None
        self.chunks = None
        self.ids_kept = None
        self.annotations = []
        self.hitsort = None
        if filename:
            logger.debug("Creating database in file")
            if os.path.isfile(filename) and (new or source):
                # remove old database if exists
                os.remove(filename)
            # sqlite database in file
            self.conn = sqlite3.connect(filename)
        else:
            # sqlite database in memory
            logger.debug("Creating database in memory")
            self.conn = sqlite3.connect(":memory:")
        c = self.conn.cursor()
        c.execute("create table sequences (name text, sequence text)")

        if source:
            self._read_from_fasta()
            c.execute("CREATE TABLE prefix_codes (prefix TEXT, N INTEGER)")
            self.conn.executemany("INSERT INTO prefix_codes (prefix, N) values (?,?)", self.prefix_codes.items())
            self.conn.commit()
        self._commit_immediately = True
        self.makeblastdb(legacy=True)
        # store some attribudes in database
        self._update_database()
        logger.info("sequences loaded")
        logger.info(self)

    def _update_database(self):
        # pass all atributes - which are either float, ind and str to extra table
        c = self.conn.cursor()
        stored_attributes = [
            ("sample_size", "integer"),
            ("_length", "integer"),
            ("_original_length", "integer"),
            ("paired", "integer"),
        ]
        columns = ", ".join(i[0] + " " + i[1] for i in stored_attributes)
        try:
            c.execute((
                "CREATE TABLE seqinfo ( {} )".format(columns)
            ))
        except sqlite3.OperationalError:
            pass # table already exists
        # get current values
        values = tuple(getattr(self, i[0]) for i in stored_attributes)
        placeholder = ", ".join(["?"] * len(values))
        c.execute("DELETE FROM seqinfo")
        c.execute("INSERT INTO seqinfo VALUES ({})".format(placeholder), values)


    def _read_from_fasta(self):
        # changed will be commited after whole file is loaded - this is faster
        self._commit_immediately = False
        if self.fasta:
            f = open(self.fasta, mode='w')
        # define sampling
        if self.sample_size:
            logger.debug(
                'sampling sequences: {} sample size'.format(self.sample_size))
            N = self.fasta_length(self.source)
            if self.paired:
                s = itertools.chain(
                    [[i, i + 1] for i in sorted(random.sample(range(1, N, 2), int(self.sample_size / 2)))])

                sample = itertools.chain(
                    [item for sublist in s for item in sublist])

            else:
                sample = itertools.chain(
                    sorted(random.sample(range(1, N + 1), self.sample_size)))
                logger.debug("sampling unpaired reads")
        else:
            # no sampling - use all
            sample = itertools.count(1)
        
        # define counter:
        if self.rename:
            if self.paired:
                counter = itertools.count(1, 0.5)
            else:
                counter = itertools.count(1, 1)
        else:
            counter = itertools.cycle([""])
        # define pairs:
        if self.paired:
            pair = itertools.cycle(['f', 'r'])
        else:
            pair = itertools.cycle([''])
        position = 0

        seq2write = next(sample)
        
        for i in Sequence.read_fasta(self.source):
            position += 1
            if position == seq2write:

                # create name:
                prefix = i.name[0:self.prefix_length] # could be empty ""
                if self.rename:
                    i.name = "{0}{1}{2}".format(
                        prefix,
                        int(next(counter)),
                        next(pair)
                    )
                if self.prefix_length:
                    if prefix in self.prefix_codes:
                        self.prefix_codes[prefix] += 1
                    else:
                        self.prefix_codes[prefix] = 1
                self[i.name] = i.seq
                if self.fasta:
                    i.write2fasta(f)
                try:
                    seq2write = next(sample)
                except StopIteration:
                    # no more sequences to sample
                    break
     

        self.conn.commit()   # save changes
        if self.fasta:
            f.close()
        self._commit_immediately = True

    def __getitem__(self, item):
        # item cannot be empty!!!
        c = self.conn.cursor()
        c.execute("select * from sequences where name= '{}'".format(item))
        # try:
        name, sequence = c.fetchone()  # should be only one
        # except TypeError:
        #    return None
        return sequence

    def __setitem__(self, item, value):
        c = self.conn.cursor()
        c.execute(
            "insert into sequences values ( '{0}', '{1}')".format(item, value))
        if self._commit_immediately:
            self.conn.commit()   # save changes

#    def __iter__(self):
#       self.c.execute("select name from sequences")
#       for i in self.c:
#            yield i[0]

    def __iter__(self):
        c = self.conn.cursor()
        c.execute("select name from sequences")
        for i in c:
            yield i[0]


#    def __iter__(self):
#        for i in range(1,len(self)):
#            yield i
    def items(self):
        c = self.conn.cursor()
        c.execute("select name, sequence from sequences")
        for name, seq in c:
            yield Sequence(seq, name)

    def toDict(self):
        c = self.conn.cursor()
        c.execute("select name, sequence from sequences")
        d = {}
        for name, seq in c:
            d[name]=seq
        return(d)

    def __len__(self):
        c = self.conn.cursor()
        c.execute("select count(*) from sequences")
        length = c.fetchone()[0]
        return(length)

    def __str__(self):
        out = ""
        c = self.conn.cursor()
        c.execute("select * from sequences limit {0}".format(MAX_PRINT))
        for n, s in c:
            out = out + "{0} : {1}\n".format(n, s)
        out = out + "...\n"
        out = out + str(len(self)) + " sequence total\n"
        return out

    def insert(self, sequence, commit=True):
        self._commit_immediately = commit
        self[sequence.name] = sequence.seq
        self._commit_immediately = True

    def keys(self):
        '''
        this will get names of sequences
        '''
        c = self.conn.cursor()
        c.execute('select name from sequences order by rowid')
        names = []
        for i in c.fetchall():
            names.append(i[0])
        return names

    def sample(self, size, seed=123):
        '''
        generator - reproducible sampling sequences
        '''
        max_index = len(self)
        # sample = seqtools.SequenceSet(source=info.input_sequences, )
        if self.paired:
            size = int(size / 2)
            step = 2
        else:
            step = 1
        random.seed(seed)
        c = self.conn.cursor()
        for i in sorted(random.sample(range(1, max_index, step), size)):
            c.execute(
                "select name, sequence from sequences where rowid={}".format(i))
            name, sequence = c.fetchone()
            yield Sequence(sequence, name)
            if self.paired:
                c.execute(
                    "select name, sequence from sequences where rowid={}".format(i + 1))
                name, sequence = c.fetchone()
                yield Sequence(sequence, name)

    def sample2db(self, db_file=None, fasta_file_name=None, size=100, seed=123):
        if (not db_file) and (not fasta_file_name):
            raise IOError("no db_file nor fasta_file_name were defined")
        # open files
        if db_file:
            db = SequenceSet(filename=db_file, source=None, new=True)
        if fasta_file_name:
            f = open(fasta_file_name, 'w')
        # write in files
        for i in self.sample(size, seed):
            if db_file:
                db.insert(i, commit=False)
            if fasta_file_name:
                i.write2fasta(f)
        # close files
        if fasta_file_name:
            f.close()
        if db_file:
            db.conn.commit()
            return db

    def save2fasta(self, fasta_file_name, keep=True, subset=None):
        if subset:
            with open(fasta_file_name, 'w') as f:
                c = self.conn.cursor()
                c.execute("select name, sequence from sequences where name in ('{}')".format(
                    "','".join(subset)))
                for name, sequence in c:
                    s = Sequence(sequence, name)
                    s.write2fasta(f)
        else:
            with open(fasta_file_name, 'w') as f:
                for i in self.items():
                    i.write2fasta(f)
            if keep:
                self.fasta = fasta_file_name

    def save_annotation(self, annotation_name, subset, dir):
        annotation_file = "{}/{}_annotation.csv".format(dir, annotation_name)
        with open(annotation_file, "w") as f:
            c = self.conn.cursor()
            c.execute(
                "select * from {} where name in ('{}')".format(
                    annotation_name, "','".join(subset)
                )
            )
            header = [description[0] for description in c.description]
            f.write("\t".join(str(j) for j in header) + "\n")
            for i in c:
                f.write("\t".join(str(j) for j in i) + "\n")
        return annotation_file

    def make_chunks(self, file_name=None, chunk_size=5000):
        '''
        split SequneceSet to chunks and save to files,
        return list of files names
        '''
        logger.debug("creating chunks from fasta file: ")
        if not file_name and self.fasta:
            file_name = self.fasta
        else:
            raise Exception('file name for chunks is not defined!')

        seqs = iter(self.items())
        fn_chunk_names = []
        for i in itertools.count():
            try:
                fn = "{0}.{1}".format(file_name, i)
                logger.info("saving chunk {}".format(fn))
                for n in range(chunk_size):
                    s = next(seqs)
                    if not n:
                        f = open(fn, 'w')
                    s.write2fasta(f)
                f.close()
                fn_chunk_names.append(fn)
            except StopIteration:
                if not f.closed:
                    f.close()
                    fn_chunk_names.append(fn)
                break
        self.chunks = fn_chunk_names
        logger.debug(self.chunks)
        return(fn_chunk_names)

    def makeblastdb(self, blastdb=None, dbtype='nucl', legacy=False, lastdb = False):
        '''
        dbtype eithe 'nucl' of 'prot'
        '''
        logger.info("creating blast database")
        # check if fasta file on disk is present
        if not self.fasta:
            try:
                self.save2fasta(blastdb)
            except TypeError:
                print("\nEither blastdb or fasta must be ",
                      "defined when making blast database!!!\n")
                raise
        else:
            blastdb = blastdb if blastdb else self.fasta

        cmd = "makeblastdb -in {0} -out {1} -dbtype {2}".format(self.fasta,
                                                                blastdb,
                                                                dbtype)
        args = shlex.split(cmd)
        if 0 == subprocess.check_call(args, stderr=sys.stdout):
            self.blastdb = blastdb
        else:
            Warning("makeblastdb failed")

        if legacy:
            logger.info("creating legacy blast database")
            prot = {'nucl': 'F', 'prot': 'T'}[dbtype]
            cmd = "formatdb -i {0} -p {1} -n {2}.legacy".format(self.fasta,
                                                                prot,
                                                                blastdb
                                                                )
            args = shlex.split(cmd)
            if subprocess.check_call(args) == 0:
                self.blastdb_legacy = "{}.legacy".format(blastdb)
            else:
                Warning("formatdb failed")
        if lastdb:
            logger.info("creating lastdb database")
            cmd = "lastdb -u NEAR  {fasta} {fasta}".format(fasta=self.fasta)
            args = shlex.split(cmd)
            if subprocess.check_call(args) == 0:
                self.lastdb = self.fasta
            else:
                Warning("formatdb failed")


    def create_hitsort(self, options, output=None):
        '''
        query is name of files to search against blastdb of SequenceSet object
        query could be multiple files
        this is inteded to be used for tabular output only
        output file must be specified,
        if more than one query is used, multiple files are created
        worker - function to execute blast search
        '''

        logger.info('running all to all blast')
        query = self.chunks
        if not output:
            output = "{}.blast".format(self.fasta)
        database = getattr(self, options.database)
        if len(query) > 1:
            output_parts = [output + "." + str(i) for i in range(len(query))]
        else:
            output_parts = [output]
            query = [query]

        parallel(_hitsort_worker, query, [database],
                 output_parts, [options])
        logger.info('all to all blast finished')
        logger.info('removing duplicates from all to all blast results')
        self._clean_hitsort(
            hitsort_files=output_parts, filtered_hitsort=output)
        ## deleting is now done during hitsort cleaning
        ##for f in output_parts:
        ##     os.unlink(f)
        self.hitsort = output
        return output

    def _clean_hitsort(self, hitsort_files, filtered_hitsort):
        ''' alterantive hitsort duplicate removal  '''

        ids = self.keys()
        Nfiles = 5 + int(sum([os.path.getsize(i)
                              for i in hitsort_files]) / 5e9)
        hitsort_parts = [
            "{0}.{1}.parts".format(filtered_hitsort, i) for i in range(Nfiles)]
        # open files
        fs = []
        for i in hitsort_parts:
            fs.append(open(i, 'w'))

        ids_destination = {}
        fs_iter = itertools.cycle(fs)
        for i in ids:
            #ids_destination[i] = fs_iter.next()
            ids_destination[i] = next(fs_iter)
        lines_out = lines_total = 0
        temp_open = False

        for i in hitsort_files:
            f = open(i, 'r')
            while True:
                line = f.readline()
                if not line:
                    break
                lines_total += 1
                line_items = line.split()
#                ids_destination[line_items[0]].write(line)
                # do not assume that it is sorted!
                ids_destination[min(line_items[0:2])].write(line)
            os.unlink(i)
        # close all files
        for i in fs:
            i.close()

        # load one by one and exclude duplicates:
        hit_out = open(filtered_hitsort, 'w')
        for i in hitsort_parts:
            f = open(i, 'r')
            hitsort_shelve = OrderedDict()
            while True:
                line = f.readline()
                if not line:
                    break
                line_items = line.split()
                key = "\t".join(sorted(line_items[0:2]))
                value = line_items
                bitscore = float(line_items[2])
                if key in hitsort_shelve:
                    if float(hitsort_shelve[key][2]) > bitscore:
                        hitsort_shelve[key] = value
                    hit_out.write("\t".join(hitsort_shelve[key]) + "\n")
                    del hitsort_shelve[key]
                    lines_out += 1
                else:
                    hitsort_shelve[key] = value
            f.close()
            for key in hitsort_shelve:
                hit_out.write("\t".join(hitsort_shelve[key]) + "\n")
                lines_out += 1
        # clean up
        for i in hitsort_parts:
            os.unlink(i)
        hit_out.close()
        return lines_total, lines_out

    def _check_database(self, database):
        '''
        database is path to fastafile, database with same prefix should exist
        it creates database if does not exist as temporal file!
        '''
        pass
        # TODO

    def remove_sequences_using_filter(self, fasta_file_filter, keep_proportion,
                                      omitted_sequences_file,
                                      kept_sequences_file):
        '''
        Run blast with fasta_file_filter against legaci blastdb of SequenceSet to detect sequences
        which are to be removed from SequencesSet.  Complete pairsare removed if paired sequences
        are used.
        '''
        ids = blast_with_filter(fasta_file_filter, self.blastdb_legacy)
        # check against database and remove sequences

        c = self.conn.cursor()
        if self.paired:
            c.execute("SELECT rowid FROM sequences WHERE name IN {}".format(
                format_query(ids)
            ))
            odd_rowids = set()
            for i in c.fetchall():
                if i[0] % 2 == 0:
                    odd_rowids.add(i[0] - 1)
                else:
                    odd_rowids.add(i[0])
            odd_rowids_keep = random.sample(odd_rowids, round(keep_proportion * len(odd_rowids)))
            c.execute("SELECT name FROM sequences WHERE rowid IN {} ORDER BY rowid".format(
                format_query(odd_rowids_keep + [i+1 for i in odd_rowids_keep])
            ))
            self.ids_kept = [i[0] for i in c.fetchall()]
            odd_rowids_omit = list(odd_rowids.difference(odd_rowids_keep))
            c.execute("SELECT name FROM sequences WHERE rowid IN {} ORDER BY rowid".format(
                format_query(odd_rowids_omit + [i+1 for i in odd_rowids_omit])
            ))
            ids_omit = [i[0] for i in c.fetchall()]

        else:
            self.ids_kept = random.sample(ids, round(keep_proportion * len(ids)))
            ids_omit = list(ids.difference(self.ids_kept))
        self.save2fasta(omitted_sequences_file, keep=False, subset=ids_omit)
        self.save2fasta(kept_sequences_file, keep=False, subset=self.ids_kept)
        self.drop_sequences(ids_omit)
        # save info about omited sequences
        c.execute("CREATE TABLE ids_kept (name TEXT)")
        c.executemany("INSERT INTO ids_kept (name) values (?)", [(i,) for i in self.ids_kept])
        return len(ids_omit)

    def drop_sequences(self, ids_to_drop):
        '''
        remove sequences from sequence set based on the ID
        '''
        self._original_length = len(self)
        c = self.conn.cursor()
        idslist = format_query(ids_to_drop)
        c.execute("CREATE TABLE omited_sequences AS SELECT * FROM sequences  WHERE name IN {}".format(idslist))
        c.execute("DELETE FROM sequences WHERE name IN {}".format(
            idslist
        ))
        # new databases must be created - with ommited sequences!
        self.save2fasta(fasta_file_name=self.fasta) ## this will replace origninal file!
        self._length = len(self)
        self.makeblastdb(legacy=True)
        self._update_database()

    def annotate(self, database, annotation_name, directory="", params=""):
        '''
        database : path to blast formated database
        method: type of search to use for annotation. it must be
        'blastn' of 'blastx'
        Annotation is save in directory also into the table with name stored in
        annotation_name variable
        '''
        logger.info("annotating reads  with {} ".format(annotation_name))
        self._check_database(database)
        if "parallelize" in params:
            parallelize = params['parallelize']
        else:
            parallelize = True
        if parallelize:
            if not self.chunks:
                self.make_chunks()
            output = [
                "{}/{}_{}".format(directory, annotation_name,os.path.basename(i)) for i in self.chunks]
            parallel(blast_worker, self.chunks, [database], output, [params])
        else:
            # no explicit parallelization using parallel
            single_output =  "{}/{}_hits".format(directory, annotation_name)
            params['args'] = params['args'].format(max_proc=get_max_proc())
            blast_worker(self.fasta, database, single_output,params)
            output = [single_output]

        # put into as new attribute and sqlite table
        c = self.conn.cursor()
        c.execute("create table {} (name text, db_id text, length real, bitscore real , pid real)".format(
            annotation_name))
        for blast_results in output:
            with open(blast_results, 'r') as f:
                for i in f:
                    items = i.split()
                    types = [str, str, float, float, float, float, float]
                    qseqid, sseqid, qlen, slen, length, pident, bitscore = [
                        f(i) for i, f in zip(items, types)]
                    c.execute('insert into {} values ("{}", "{}", "{}", "{}", "{}")'.format(
                        annotation_name, qseqid, sseqid, length, bitscore, pident))
        self.conn.commit()
        self.annotations.append(annotation_name)

    @staticmethod
    def fasta_length(fasta_file_name):
        '''
        count number of sequences, 
        '''
        with open(fasta_file_name, 'r') as f:
            N = 0
            for i in f:
                if i.strip()[0] == ">":
                    N += 1
            return N


# test
if __name__ == "__main__":
    # get root directory

    MAIN_DIR = os.path.realpath(
        os.path.dirname(os.path.realpath(__file__)) + "/..")
    print("MAIN_DIR:")
    print(MAIN_DIR)
    TMP_DIR = "{}/tmp".format(MAIN_DIR)
    TEST_DATA_DIR = "{}/test_data".format(MAIN_DIR)

    # some data:
    fn = "{}/pair_reads.fasta".format(TEST_DATA_DIR)
    os.chdir(TMP_DIR)
    # s = SequenceSet(filename='sqlite.db')
    print("number of sequences in fasta file:")
    print(SequenceSet.fasta_length(fn))

    print("loading sequences...")
    s = SequenceSet(source=fn, paired=True, rename=False,
                    filename='sqlite.db', prefix_length=4, fasta='sample2.fasta')
    print(s)

    print("creating blast database")
    s.makeblastdb()

    print("creating chunks from fasta file: ")
    file_list = s.make_chunks("fasta_chunks", chunk_size=500)
    print(file_list)

    print('creating hitsort')
    s.create_hitsort(file_list, "hitsort")

    print("saving to fasta file")
    s.save2fasta('sampleX.fasta', keep=False)

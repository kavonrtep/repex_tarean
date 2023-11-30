#!/usr/bin/env python3
'''
This module is mainly for large graph (e.i hitsort) storage, parsing and for clustering
'''
import os
import sys
import sqlite3
import time
import subprocess
import logging
from collections import defaultdict
import collections
import operator
import math
import random
import itertools
import config
from lib import r2py
from lib.utils import FilePath
from lib.parallel.parallel import parallel2 as parallel
REQUIRED_VERSION = (3, 4)
MAX_BUFFER_SIZE = 100000
if sys.version_info < REQUIRED_VERSION:
    raise Exception("\n\npython 3.4 or higher is required!\n")
LOGGER = logging.getLogger(__name__)


def dfs(start, graph):
    """
    helper function for cluster merging.
    Does depth-first search, returning a set of all nodes seen.
    Takes: a graph in node --> [neighbors] form.
    """
    visited, worklist = set(), [start]

    while worklist:
        node = worklist.pop()
        if node not in visited:
            visited.add(node)
            # Add all the neighbors to the worklist.
            worklist.extend(graph[node])

    return visited


def graph_components(edges):
    """
    Given a graph as a list of edges, divide the nodes into components.
    Takes a list of pairs of nodes, where the nodes are integers.
    """

    # Construct a graph (mapping node --> [neighbors]) from the edges.
    graph = defaultdict(list)
    nodes = set()

    for v1, v2 in edges:
        nodes.add(v1)
        nodes.add(v2)

        graph[v1].append(v2)
        graph[v2].append(v1)

    # Traverse the graph to find the components.
    components = []

    # We don't care what order we see the nodes in.
    while nodes:
        component = dfs(nodes.pop(), graph)
        components.append(component)

        # Remove this component from the nodes under consideration.
        nodes -= component

    return components


class Graph():
    '''
    create Graph object stored in sqlite database, either in memory or on disk
    structure of table is:
    V1 V2 weigth12
    V2 V3 weight23
    V4 V5 weight45
    ...
    ...
    !! this is undirected simple graph - duplicated edges must
    be removed on graph creation

    '''
    # seed for random number generator - this is necessary for reproducibility between runs
    seed = '123'

    def __init__(self,
                 source=None,
                 filename=None,
                 new=False,
                 paired=True,
                 seqids=None):
        '''
        filename : fite where to store database, if not defined it is stored in memory
        source : ncol file from which describe graph
        new : if false and source is not define graph can be loaded from database (filename)

        vertices_name must be in correcti order!!!
        '''

        self.filename = filename
        self.source = source
        self.paired = paired
        # path to indexed graph - will be set later
        self.indexed_file = None
        self._cluster_list = None
        # these two attributes are set after clustering
        # communities before merging
        self.graph_2_community0 = None
        # communities after merging
        self.graph_2_community = None
        self.number_of_clusters = None
        self.binary_file = None
        self.cluster_sizes = None
        self.graph_tree = None
        self.graph_tree_log = None
        self.weights_file = None

        if filename:
            if os.path.isfile(filename) and (new or source):
                os.remove(filename)
            self.conn = sqlite3.connect(filename)
        else:
            self.conn = sqlite3.connect(":memory:")
        c = self.conn.cursor()

        c.execute("PRAGMA page_size=8192")
        c.execute("PRAGMA cache_size = 2000000 ")  # this helps

        try:
            c.execute((
                "create table graph (v1 integer, v2 integer, weight integer, "
                "pair integer, v1length integer, v1start  integer, v1end  integer, "
                "v2length  integer, v2start  integer, v2end  integer, pid  integer,"
                "evalue real, strand text )"))
        except sqlite3.OperationalError:
            pass  # table already exist
        else:
            c.execute(
                "create table vertices (vertexname text primary key, vertexindex integer)")
        tables = sorted(c.execute(
            "SELECT name FROM sqlite_master WHERE type='table'").fetchall())

        if not [('graph', ), ('vertices', )] == tables:
            raise Exception("tables for sqlite for graph are not correct")

        if source:
            self._read_from_hitsort()

        if paired and seqids:
            # vertices must be defined - create graph of paired reads:
            # last character must disinguish pair
            c.execute((
                "create table pairs (basename, vertexname1, vertexname2,"
                "v1 integer, v2 integer, cluster1 integer, cluster2 integer)"))
            buffer = []
            for i, k in zip(seqids[0::2], seqids[1::2]):
                assert i[:-1] == k[:-1], "problem with pair reads ids"
                # some vertices are not in graph - singletons
                try:
                    index1 = self.vertices[i]
                except KeyError:
                    index1 = -1

                try:
                    index2 = self.vertices[k]
                except KeyError:
                    index2 = -1

                buffer.append((i[:-1], i, k, index1, index2))

            self.conn.executemany(
                "insert into pairs (basename, vertexname1, vertexname2, v1, v2) values (?,?,?,?,?)",
                buffer)
            self.conn.commit()

    def _read_from_hitsort(self):

        c = self.conn.cursor()
        c.execute("delete from graph")
        buffer = []
        vertices = {}
        counter = 0
        v_count = 0
        with open(self.source, 'r') as f:
            for i in f:
                edge_index = {}
                items = i.split()
                # get or insert vertex index
                for vn in items[0:2]:
                    if vn not in vertices:
                        vertices[vn] = v_count
                        edge_index[vn] = v_count
                        v_count += 1
                    else:
                        edge_index[vn] = vertices[vn]
                if self.paired:
                    pair = int(items[0][:-1] == items[1][:-1])
                else:
                    pair = 0
                buffer.append(((edge_index[items[0]], edge_index[items[1]],
                                items[2], pair) + tuple(items[3:])))
                if len(buffer) == MAX_BUFFER_SIZE:
                    counter += 1
                    self.conn.executemany(
                        "insert or ignore into graph values (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        buffer)
                    buffer = []
        if buffer:
            self.conn.executemany(
                "insert or ignore into graph values (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                buffer)

        self.conn.commit()
        self.vertices = vertices
        self.vertexid2name = {
            vertex: index
            for index, vertex in vertices.items()
        }
        self.vcount = len(vertices)
        c = self.conn.cursor()
        c.execute("select count(*) from graph")
        self.ecount = c.fetchone()[0]
        # fill table of vertices
        self.conn.executemany("insert into vertices values (?,?)",
                              vertices.items())
        self.conn.commit()

    def save_indexed_graph(self, file=None):
        if not file:
            self.indexed_file = "{}.int".format(self.source)
        else:
            self.indexed_file = file
        c = self.conn.cursor()
        with open(self.indexed_file, 'w') as f:
            out = c.execute('select v1,v2,weight from graph')
            for v1, v2, weight in out:
                f.write('{}\t{}\t{}\n'.format(v1, v2, weight))

    def get_subgraph(self, vertices):
        pass

    def _levels(self):
        with open(self.graph_tree_log, 'r') as f:
            levels = -1
            for i in f:
                if i[:5] == 'level':
                    levels += 1
        return levels

    def _reindex_community(self, id2com):
        '''
        reindex community and superclusters so that biggest cluster is no.1
        '''
        self.conn.commit()
        _, community, supercluster = zip(*id2com)
        (cluster_index, frq, self.cluster_sizes,
         self.number_of_clusters) = self._get_index_and_frequency(community)

        supercluster_index, sc_frq, _, _ = self._get_index_and_frequency(
            supercluster)
        id2com_reindexed = []

        for i, _ in enumerate(id2com):
            id2com_reindexed.append((id2com[i][0], id2com[i][1], frq[
                i], cluster_index[i], supercluster_index[i], sc_frq[i]))
        return id2com_reindexed

    @staticmethod
    def _get_index_and_frequency(membership):
        frequency_table = collections.Counter(membership)
        frequency_table_sorted = sorted(frequency_table.items(),
                                        key=operator.itemgetter(1),
                                        reverse=True)
        frq = []
        for i in membership:
            frq.append(frequency_table[i])
        rank = {}
        index = 0
        for comm, _ in frequency_table_sorted:
            index += 1
            rank[comm] = index
        cluster_index = [rank[i] for i in membership]
        cluster_sizes = [i[1] for i in frequency_table_sorted]
        number_of_clusters = len(frequency_table)
        return [cluster_index, frq, cluster_sizes, number_of_clusters]

    def louvain_clustering(self, merge_threshold=0, cleanup=False):
        '''
        input - graph
        output - list of clusters
        executables path ??
        '''
        LOGGER.info("converting hitsort to binary format")
        self.binary_file = "{}.bin".format(self.indexed_file)
        self.weights_file = "{}.weight".format(self.indexed_file)
        self.graph_tree = "{}.graph_tree".format(self.indexed_file)
        self.graph_tree_log = "{}.graph_tree_log".format(self.indexed_file)
        self.graph_2_community0 = "{}.graph_2_community0".format(
            self.indexed_file)
        self._cluster_list = None
        self.graph_2_community = "{}.graph_2_community".format(
            self.indexed_file)
        print(["louvain_convert", "-i", self.indexed_file, "-o",
               self.binary_file, "-w", self.weights_file])
        subprocess.check_call(
            ["louvain_convert", "-i", self.indexed_file, "-o",
             self.binary_file, "-w", self.weights_file],
            timeout=None)

        gt = open(self.graph_tree, 'w')
        gtl = open(self.graph_tree_log, 'w')
        LOGGER.info("running louvain clustering...")
        subprocess.check_call(
            ["louvain_community", self.binary_file, "-l", "-1", "-w",
             self.weights_file, "-v ", "-s", self.seed],
            stdout=gt,
            stderr=gtl,
            timeout=None)
        gt.close()
        gtl.close()

        LOGGER.info("creating list of cummunities")
        gt2c = open(self.graph_2_community0, 'w')
        subprocess.check_call(
            ['louvain_hierarchy', self.graph_tree, "-l", str(self._levels())],
            stdout=gt2c)
        gt2c.close()
        if merge_threshold and self.paired:
            com2newcom = self.find_superclusters(merge_threshold)
        elif self.paired:
            com2newcom = self.find_superclusters(config.SUPERCLUSTER_THRESHOLD)
        else:
            com2newcom = {}
        # merging of clusters, creatting superclusters
        LOGGER.info("mergings clusters based on mate-pairs ")
        # modify self.graph_2_community file
        # rewrite graph2community
        with open(self.graph_2_community0, 'r') as fin:
            with open(self.graph_2_community, 'w') as fout:
                for i in fin:
                    # write  graph 2 community file in format:
                    # id communityid supeclusterid
                    # if merging - community and superclustwers are identical
                    vi, com = i.split()
                    if merge_threshold:
                        ## mergin
                        if int(com) in com2newcom:
                            fout.write("{} {} {}\n".format(vi, com2newcom[int(
                                com)], com2newcom[int(com)]))
                        else:
                            fout.write("{} {} {}\n".format(vi, com, com))
                    else:
                        ## superclusters
                        if int(com) in com2newcom:
                            fout.write("{} {} {}\n".format(vi, com, com2newcom[
                                int(com)]))
                        else:
                            fout.write("{} {} {}\n".format(vi, com, com))

        LOGGER.info("loading communities into database")
        c = self.conn.cursor()
        c.execute(("create table communities (vertexindex integer primary key,"
                   "community integer, size integer, cluster integer, "
                   "supercluster integer, supercluster_size integer)"))
        id2com = []
        with open(self.graph_2_community, 'r') as f:
            for i in f:
                name, com, supercluster = i.split()
                id2com.append((name, com, supercluster))
        id2com_reindexed = self._reindex_community(id2com)
        c.executemany("insert into communities values (?,?,?,?,?,?)",
                      id2com_reindexed)
        #create table of superclusters - clusters
        c.execute(("create table superclusters as "
                   "select distinct supercluster, supercluster_size, "
                   "cluster, size from communities;"))
        # create view id-index-cluster
        c.execute(
            ("CREATE VIEW vertex_cluster AS SELECT vertices.vertexname,"
             "vertices.vertexindex, communities.cluster, communities.size"
             " FROM vertices JOIN communities USING (vertexindex)"))
        self.conn.commit()

        # add clustering infor to graph
        LOGGER.info("updating graph table")
        t0 = time.time()

        c.execute("alter table graph add c1 integer")
        c.execute("alter table graph add c2 integer")
        c.execute(("update graph set c1 = (select cluster  FROM communities "
                   "where communities.vertexindex=graph.v1)"))
        c.execute(
            ("update graph set c2 = (select cluster  FROM communities where "
             "communities.vertexindex=graph.v2)"))
        self.conn.commit()
        t1 = time.time()
        LOGGER.info("updating graph table - done in {} seconds".format(t1 -
                                                                       t0))

        # identify similarity  connections between clusters
        c.execute(
            "create table cluster_connections as SELECT c1,c2 , count(*) FROM (SELECT c1, c2  FROM graph WHERE c1>c2 UNION ALL SELECT c2 as c1, c1  as c2 FROM graph WHERE c2>c1) GROUP BY c1, c2")
        # TODO - remove directionality - summarize -

        # add cluster identity to pairs table

        if self.paired:
            LOGGER.info("analyzing pairs ")
            t0 = time.time()
            c.execute(
                "UPDATE pairs SET cluster1=(SELECT cluster FROM communities WHERE communities.vertexindex=pairs.v1)")
            t1 = time.time()
            LOGGER.info(
                "updating pairs table - cluster1 - done in {} seconds".format(
                    t1 - t0))

            t0 = time.time()
            c.execute(
                "UPDATE pairs SET cluster2=(SELECT cluster FROM communities WHERE communities.vertexindex=pairs.v2)")
            t1 = time.time()
            LOGGER.info(
                "updating pairs table - cluster2 - done in {} seconds".format(
                    t1 - t0))
            # reorder records

            t0 = time.time()
            c.execute(
                "UPDATE pairs SET cluster1=cluster2, cluster2=cluster1, vertexname1=vertexname2,vertexname2=vertexname1 where cluster1<cluster2")
            t1 = time.time()
            LOGGER.info("sorting - done in {} seconds".format(t1 - t0))

            t0 = time.time()
            c.execute(
                "create table cluster_mate_connections as select cluster1 as c1, cluster2 as c2, count(*) as N, group_concat(basename) as ids from pairs where cluster1!=cluster2 group by cluster1, cluster2;")
            t1 = time.time()
            LOGGER.info(
                "creating cluster_mate_connections table - done in {} seconds".format(
                    t1 - t0))
            # summarize
            t0 = time.time()
            self._calculate_pair_bond()
            t1 = time.time()
            LOGGER.info(
                "calculating cluster pair bond - done in {} seconds".format(
                    t1 - t0))
            t0 = time.time()
        else:
            # not paired - create empty tables
            self._add_empty_tables()

        self.conn.commit()
        t1 = time.time()
        LOGGER.info("commiting changes - done in {} seconds".format(t1 - t0))

        if cleanup:
            LOGGER.info("cleaning clustering temp files")
            os.unlink(self.binary_file)
            os.unlink(self.weights_file)
            os.unlink(self.graph_tree)
            os.unlink(self.graph_tree_log)
            os.unlink(self.graph_2_community0)
            os.unlink(self.graph_2_community)
            os.unlink(self.indexed_file)
            self.binary_file = None
            self.weights_file = None
            self.graph_tree = None
            self.graph_tree_log = None
            self.graph_2_community0 = None
            self.graph_2_community = None
            self.indexed_file = None

        # calcultate k

    def find_superclusters(self, merge_threshold):
        '''Find superclusters from clustering based on paired reads '''
        clsdict = {}
        with open(self.graph_2_community0, 'r') as f:
            for i in f:
                vi, com = i.split()
                if com in clsdict:
                    clsdict[com] += [self.vertexid2name[int(vi)][0:-1]]
                else:
                    clsdict[com] = [self.vertexid2name[int(vi)][0:-1]]
        # remove all small clusters - these will not be merged:
        small_cls = []
        for i in clsdict:
            if len(clsdict[i]) < config.MINIMUM_NUMBER_OF_READS_FOR_MERGING:
                small_cls.append(i)
        for i in small_cls:
            del clsdict[i]
        pairs = []
        for i, j in itertools.combinations(clsdict, 2):
            s1 = set(clsdict[i])
            s2 = set(clsdict[j])
            wgh = len(s1 & s2)
            if wgh < config.MINIMUM_NUMBER_OF_SHARED_PAIRS_FOR_MERGING:
                continue
            else:
                n1 = len(s1) * 2 - len(clsdict[i])
                n2 = len(s2) * 2 - len(clsdict[j])
                k = 2 * wgh / (n1 + n2)
            if k > merge_threshold:
                pairs.append((int(i), int(j)))
        # find connected commponents - will be merged
        cls2merge = graph_components(pairs)
        com2newcom = {}
        for i in cls2merge:
            newcom = min(i)
            for j in i:
                com2newcom[j] = newcom
        return com2newcom

    def adjust_cluster_size(self, proportion_kept, ids_kept):
        LOGGER.info("adjusting cluster sizes")
        c = self.conn.cursor()
        c.execute("ALTER TABLE superclusters ADD COLUMN size_uncorrected INTEGER")
        c.execute("UPDATE superclusters SET size_uncorrected=size")
        if ids_kept:
            ids_kept_set = set(ids_kept)
            ratio = (1 - proportion_kept)/proportion_kept
            for cl, size in c.execute("SELECT cluster,size FROM superclusters"):
                ids = self.get_cluster_reads(cl)
                ovl_size = len(ids_kept_set.intersection(ids))
                size_adjusted = int(len(ids) + ovl_size * ratio)
                if size_adjusted > size:
                    c.execute("UPDATE superclusters SET size=? WHERE cluster=?",
                              (size_adjusted, cl))
        self.conn.commit()
        LOGGER.info("adjusting cluster sizes - done")

    def export_cls(self, path):
        with open(path, 'w') as f:
            for i in range(1, self.number_of_clusters + 1):
                ids = self.get_cluster_reads(i)
                f.write(">CL{}\t{}\n".format(i, len(ids)))
                f.write("\t".join(ids))
                f.write("\n")

    def _calculate_pair_bond(self):
        c = self.conn.cursor()
        out = c.execute("select c1, c2, ids from cluster_mate_connections")
        buffer = []
        for c1, c2, ids in out:
            w = len(set(ids.split(",")))
            n1 = len(set([i[:-1] for i in self.get_cluster_reads(c1)
                          ])) * 2 - len(self.get_cluster_reads(c1))
            n2 = len(set([i[:-1] for i in self.get_cluster_reads(c2)
                          ])) * 2 - len(self.get_cluster_reads(c2))
            buffer.append((c1, c2, n1, n2, w, 2 * w / (n1 + n2)))
        c.execute(
            "CREATE TABLE cluster_mate_bond (c1 INTEGER, c2 INTEGER, n1 INTEGER, n2 INTEGER, w INTEGER, k FLOAT)")
        c.executemany(" INSERT INTO cluster_mate_bond values (?,?,?,?,?,?)",
                      buffer)

    def _add_empty_tables(self):
        '''This is used with reads that are not paired
        - it creates empty mate tables, this is necessary for
        subsequent reporting to work corectly '''
        c = self.conn.cursor()
        c.execute(("CREATE TABLE cluster_mate_bond (c1 INTEGER, c2 INTEGER, "
                   "n1 INTEGER, n2 INTEGER, w INTEGER, k FLOAT)"))
        c.execute(
            "CREATE TABLE cluster_mate_connections (c1 INTEGER, c2 INTEGER, N INTEGER, ids TEXT) ")

    def get_cluster_supercluster(self, cluster):
        '''Get supercluster id for suplied cluster '''
        c = self.conn.cursor()
        out = c.execute(
            'SELECT supercluster FROM communities WHERE cluster="{0}" LIMIT 1'.format(
                cluster))
        sc = out.fetchone()[0]
        return sc

    def get_cluster_reads(self, cluster):

        if self._cluster_list:
            return self._cluster_list[str(cluster)]
        else:
            # if queried first time
            c = self.conn.cursor()
            out = c.execute("select cluster, vertexname from vertex_cluster")
            cluster_list = collections.defaultdict(list)
            for clusterindex, vertexname in out:
                cluster_list[str(clusterindex)].append(vertexname)
            self._cluster_list = cluster_list
            return self._cluster_list[str(cluster)]

        
    def extract_cluster_blast(self, path, index, ids=None):
        ''' Extract blast for cluster and save it to path
        return number of blast lines ( i.e. number of graph edges E)
        if ids is specified , only subset of blast is used'''
        c = self.conn.cursor()
        if ids:
            vertexindex = (
                "select vertexindex from vertices "
                "where vertexname in ({})").format('"' + '","'.join(ids) + '"')

            out = c.execute(("select * from graph where c1={0} and c2={0}"
                             " and v1 in ({1}) and v2 in ({1})").format(
                                 index, vertexindex))
        else:
            out = c.execute(
                "select * from graph where c1={0} and c2={0}".format(index))
        E = 0
        N = len(self.get_cluster_reads(index))
        with open(path, 'w') as f:
            for i in out:
                print(self.vertexid2name[i[0]],
                      self.vertexid2name[
                          i[1]],
                      i[2],
                      *i[4:13],
                      sep='\t',
                      file=f)
                E += 1
        return E

    def export_clusters_files_multiple(self,
                                       min_size,
                                       directory,
                                       sequences=None,
                                       tRNA_database_path=None,
                                       satellite_model_path=None):
        def load_fun(N, E):
            ''' estimate mem usage from graph size and density'''
            NE = math.log(float(N) * float(E), 10)
            if NE > 11.5:
                return 1
            if NE > 11:
                return 0.9
            if NE > 10:
                return 0.4
            if NE > 9:
                return 0.2
            if NE > 8:
                return 0.07
            return 0.02

        def estimate_sample_size(NV, NE, maxv, maxe):
            ''' estimat suitable sampling based on the graph density
            NV,NE is |V| and |E| of the graph
            maxv, maxe are maximal |V| and |E|'''

            d = (2 * NE) / (NV * (NV - 1))
            eEst = (maxv * (maxv - 1) * d) / 2
            nEst = (d + math.sqrt(d**2 + 8 * d * maxe)) / (2 * d)
            if eEst >= maxe:
                N = int(nEst)
            if nEst >= maxv:
                N = int(maxv)
            return N

        clusterindex = 1
        cluster_input_args = []
        ppn = []
        # is is comparative analysis?
        if sequences.prefix_length:
            self.conn.execute("CREATE TABLE comparative_counts (clusterindex INTEGER,"
                    + ", ".join(["[{}] INTEGER".format(i) for i in sequences.prefix_codes.keys()]) + ")")
                    # do for comparative analysis

            for cl in range(self.number_of_clusters):
                prefix_codes = dict((key, 0) for key in sequences.prefix_codes.keys())
                for i in self.get_cluster_reads(cl):
                    prefix_codes[i[0:sequences.prefix_length]] += 1
                header = ", ".join(["[" + str(i) + "]" for i in prefix_codes.keys()])
                values = ", ".join([str(i) for i in prefix_codes.values()])
                self.conn.execute(
                    "INSERT INTO comparative_counts (clusterindex, {}) VALUES ({}, {})".format(
                        header, cl, values))
        else:
            prefix_codes = {}

        while True:
            read_names = self.get_cluster_reads(clusterindex)
            supercluster = self.get_cluster_supercluster(clusterindex)
            N = len(read_names)
            print("sequences.ids_kept -2 ")
            print(sequences.ids_kept)
            if sequences.ids_kept:
                N_adjusted = round(len(set(sequences.ids_kept).intersection(read_names)) *
                                   ((1 - config.FILTER_PROPORTION_OF_KEPT) /
                                    config.FILTER_PROPORTION_OF_KEPT) + N)
            else:
                N_adjusted = N
            if N < min_size:
                break
            else:
                LOGGER.info("exporting cluster {}".format(clusterindex))
                blast_file = "{dir}/dir_CL{i:04}/hitsort_part.csv".format(
                    dir=directory, i=clusterindex)
                cluster_dir = "{dir}/dir_CL{i:04}".format(dir=directory,
                                                          i=clusterindex)
                fasta_file = "{dir}/reads_selection.fasta".format(dir=cluster_dir)
                fasta_file_full = "{dir}/reads.fasta".format(dir=cluster_dir)

                os.makedirs(os.path.dirname(blast_file), exist_ok=True)
                E = self.extract_cluster_blast(index=clusterindex,
                                               path=blast_file)
                # check if blast must be sampled
                n_sample = estimate_sample_size(NV=N,
                                                NE=E,
                                                maxv=config.CLUSTER_VMAX,
                                                maxe=config.CLUSTER_EMAX)
                LOGGER.info("directories created..")
                if n_sample < N:
                    LOGGER.info(("cluster is too large - sampling.."
                                 "original size: {N}\n"
                                 "sample size: {NS}\n"
                                 "").format(N=N, NS=n_sample))
                    random.seed(self.seed)
                    read_names_sample = random.sample(read_names, n_sample)
                    LOGGER.info("reads id sampled...")
                    blast_file_sample = "{dir}/dir_CL{i:04}/blast_sample.csv".format(
                        dir=directory, i=clusterindex)
                    E_sample = self.extract_cluster_blast(
                        index=clusterindex,
                        path=blast_file,
                        ids=read_names_sample)
                    LOGGER.info("numner of edges in sample: {}".format(
                        E_sample))
                    sequences.save2fasta(fasta_file, subset=read_names_sample)
                    sequences.save2fasta(fasta_file_full, subset=read_names)

                else:
                    read_names_sample = None
                    E_sample = None
                    blast_file_sample = None
                    n_sample = None
                    sequences.save2fasta(fasta_file_full, subset=read_names)
                    ## TODO - use symlink instead of :
                    sequences.save2fasta(fasta_file, subset=read_names)
                # export individual annotations tables:
                # annotation is always for full cluster
                LOGGER.info("exporting cluster annotation")
                annotations = {}
                annotations_custom = {}
                for n in sequences.annotations:
                    print("sequences.annotations:", n)
                    if n.find("custom_db") == 0:
                        print("custom")
                        annotations_custom[n] = sequences.save_annotation(
                            annotation_name=n,
                            subset=read_names,
                            dir=cluster_dir)
                    else:
                        print("built in")
                        annotations[n] = sequences.save_annotation(
                            annotation_name=n,
                            subset=read_names,
                            dir=cluster_dir)

                cluster_input_args.append([
                    n_sample, N,N_adjusted, blast_file, fasta_file, fasta_file_full,
                    clusterindex, supercluster, self.paired,
                    tRNA_database_path, satellite_model_path, sequences.prefix_codes,
                    prefix_codes, annotations, annotations_custom
                ])
                clusterindex += 1
                ppn.append(load_fun(N, E))

                    

        self.conn.commit()

        # run in parallel:
        # reorder jobs based on the ppn:
        cluster_input_args = [
            x
            for (y, x) in sorted(
                zip(ppn, cluster_input_args),
                key=lambda pair: pair[0],
                reverse=True)
        ]
        ppn = sorted(ppn, reverse=True)
        LOGGER.info("creating clusters in parallel")
        clusters_info = parallel(Cluster,
                                 *[list(i) for i in zip(*cluster_input_args)],
                                 ppn=ppn)
        # sort it back:
        clusters_info = sorted(clusters_info, key=lambda cl: cl.index)
        return clusters_info


class Cluster():
    ''' store and show information about cluster properties '''

    def __init__(self,
                 size,
                 size_real,
                 size_adjusted,
                 blast_file,
                 fasta_file,
                 fasta_file_full,
                 index,
                 supercluster,
                 paired,
                 tRNA_database_path,
                 satellite_model_path,
                 all_prefix_codes,
                 prefix_codes,
                 annotations,
                 annotations_custom={},
                 loop_index_threshold=0.7,
                 pair_completeness_threshold=0.40,
                 loop_index_unpaired_threshold=0.85):
        if size:
            # cluster was scaled down
            self.size = size
            self.size_real = size_real
        else:
            self.size = self.size_real = size_real
        self.size_adjusted = size_adjusted
        self.filtered = True if size_adjusted != size_real else False
        self.all_prefix_codes = all_prefix_codes.keys
        self.prefix_codes = prefix_codes
        self.dir = FilePath(os.path.dirname(blast_file))
        self.blast_file = FilePath(blast_file)
        self.fasta_file = FilePath(fasta_file)
        self.fasta_file_full = FilePath(fasta_file_full)
        self.index = index
        self.assembly_files = {}
        self.ltr_detection = None
        self.supercluster = supercluster
        self.annotations_files = annotations
        self.annotations_files_custom = annotations_custom
        self.annotations_summary, self.annotations_table = self._summarize_annotations(
            annotations, size_real)
        # add annotation
        if len(annotations_custom):
            self.annotations_summary_custom, self.annotations_custom_table = self._summarize_annotations(
                annotations_custom, size_real)
        else:
            self.annotations_summary_custom, self.annotations_custom_table = "", ""

        self.paired = paired
        self.graph_file = FilePath("{0}/graph_layout.GL".format(self.dir))
        self.directed_graph_file = FilePath(
            "{0}/graph_layout_directed.RData".format(self.dir))
        self.fasta_oriented_file = FilePath("{0}/reads_selection_oriented.fasta".format(
            self.dir))
        self.image_file = FilePath("{0}/graph_layout.png".format(self.dir))
        self.image_file_tmb = FilePath("{0}/graph_layout_tmb.png".format(self.dir))
        self.html_report_main = FilePath("{0}/index.html".format(self.dir))
        self.html_report_files = FilePath("{0}/html_files".format(self.dir))
        self.supercluster_best_hit = "NA"
        TAREAN = r2py.R(config.RSOURCE_tarean)
        LOGGER.info("creating graph no.{}".format(self.index))
        # if FileType muast be converted to str for rfunctions
        graph_info = eval(
            TAREAN.mgblast2graph(
                self.blast_file,
                seqfile=self.fasta_file,
                seqfile_full=self.fasta_file_full,
                graph_destination=self.graph_file,
                directed_graph_destination=self.directed_graph_file,
                oriented_sequences=self.fasta_oriented_file,
                image_file=self.image_file,
                image_file_tmb=self.image_file_tmb,
                repex=True,
                paired=self.paired,
                satellite_model_path=satellite_model_path,
                maxv=config.CLUSTER_VMAX,
                maxe=config.CLUSTER_EMAX)
        )
        print(graph_info)
        self.ecount = graph_info['ecount']
        self.vcount = graph_info['vcount']
        self.loop_index = graph_info['loop_index']
        self.pair_completeness = graph_info['pair_completeness']
        self.orientation_score = graph_info['escore']
        self.satellite_probability = graph_info['satellite_probability']
        self.satellite = graph_info['satellite']
        # for paired reads:
        cond1 = (self.paired and self.loop_index > loop_index_threshold and
                 self.pair_completeness > pair_completeness_threshold)
        # no pairs
        cond2 = ((not self.paired) and
                 self.loop_index > loop_index_unpaired_threshold)
        if (cond1 or cond2) and config.ARGS.options.name != "oxford_nanopore":
            self.putative_tandem = True
            self.dir_tarean = FilePath("{}/tarean".format(self.dir))
            lock_file = self.dir + "../lock"
            out = eval(
                TAREAN.tarean(input_sequences=self.fasta_oriented_file,
                              output_dir=self.dir_tarean,
                              CPU=1,
                              reorient_reads=False,
                              tRNA_database_path=tRNA_database_path,
                              lock_file=lock_file)
            )
            self.html_tarean = FilePath(out['htmlfile'])
            self.tarean_contig_file = out['tarean_contig_file']
            self.TR_score = out['TR_score']
            self.TR_monomer_length = out['TR_monomer_length']
            self.TR_consensus = out['TR_consensus']
            self.pbs_score = out['pbs_score']
            self.max_ORF_length = out['orf_l']
            if (out['orf_l'] > config.ORF_THRESHOLD or
                    out['pbs_score'] > config.PBS_THRESHOLD):
                self.tandem_rank = 3
            elif self.satellite:
                self.tandem_rank = 1
            else:
                self.tandem_rank = 2
            # some tandems could be rDNA genes - this must be check
            # by annotation
            if self.annotations_table:
                rdna_score = 0
                contamination_score = 0
                for i in self.annotations_table:
                    if 'rDNA/' in i[0]:
                        rdna_score += i[1]
                    if 'contamination' in i[0]:
                        contamination_score += i[1]
                if rdna_score > config.RDNA_THRESHOLD:
                    self.tandem_rank = 4
                if contamination_score > config.CONTAMINATION_THRESHOLD:
                    self.tandem_rank = 0  # other

            # by custom annotation - castom annotation has preference
            if self.annotations_custom_table:
                print("custom table searching")
                rdna_score = 0
                contamination_score = 0
                print(self.annotations_custom_table)
                for i in self.annotations_custom_table:
                    if 'rDNA' in i[0]:
                        rdna_score += i[1]
                    if 'contamination' in i[0]:
                        contamination_score += i[1]
                if rdna_score > 0:
                    self.tandem_rank = 4
                if contamination_score > config.CONTAMINATION_THRESHOLD:
                    self.tandem_rank = 0  # other

        else:
            self.putative_tandem = False
            self.dir_tarean = None
            self.html_tarean = None
            self.TR_score = None
            self.TR_monomer_length = None
            self.TR_consensus = None
            self.pbs_score = None
            self.max_ORF_length = None
            self.tandem_rank = 0
            self.tarean_contig_file = None

    def __str__(self):
        out = [
            "cluster no {}:".format(self.index),
            "Number of vertices  : {}".format(self.size),
            "Number of edges     : {}".format(self.ecount),
            "Loop index          : {}".format(self.loop_index),
            "Pair completeness   : {}".format(self.pair_completeness),
            "Orientation score   : {}".format(self.orientation_score)
        ]
        return "\n".join(out)

    def listing(self, asdict=True):
        ''' convert attributes to dictionary for printing purposes'''
        out = {}
        for i in dir(self):
            # do not show private
            if i[:2] != "__":
                value = getattr(self, i)
                if not callable(value):
                    # for dictionary
                    if isinstance(value, dict):
                        for k in value:
                            out[i + "_" + k] = value[k]
                    else:
                        out[i] = value
        if asdict:
            return out
        else:
            return {'keys': list(out.keys()), 'values': list(out.values())}

    def detect_ltr(self, trna_database):
        '''detection of ltr in assembly files, output of analysis is stored in file'''
        CREATE_ANNOTATION = r2py.R(config.RSOURCE_create_annotation, verbose=False)
        if self.assembly_files['{}.{}.ace']:
            ace_file = self.assembly_files['{}.{}.ace']
            print(ace_file, "running LTR detection")
            fout = "{}/{}".format(self.dir, config.LTR_DETECTION_FILES['BASE'])
            subprocess.check_call([
                config.LTR_DETECTION,
                '-i', ace_file,
                '-o', fout,
                '-p', trna_database])
            # evaluate LTR presence
            fn = "{}/{}".format(self.dir, config.LTR_DETECTION_FILES['PBS_BLAST'])
            self.ltr_detection = CREATE_ANNOTATION.evaluate_LTR_detection(fn)


    @staticmethod
    def _summarize_annotations(annotations_files, size):
        ''' will tabulate annotation results '''
        # TODO
        summaries = {}
        # weight is in percentage
        weight = 100 / size
        for i in annotations_files:
            with open(annotations_files[i]) as f:
                header = f.readline().split()
                id_index = [
                    i for i, item in enumerate(header) if item == "db_id"
                ][0]
                for line in f:
                    classification = line.split()[id_index].split("#")[1]
                    if classification in summaries:
                        summaries[classification] += weight
                    else:
                        summaries[classification] = weight
        # format summaries for printing
        annotation_string = ""
        annotation_table = []
        for i in sorted(summaries.items(), key=lambda x: x[1], reverse=True):
            ## hits with smaller proportion are not shown!
            if i[1] > 0.1:
                if i[1] > 1:
                    annotation_string += "<b>{1:.2f}% {0}</b>\n".format(*i)
                else:
                    annotation_string += "{1:.2f}% {0}\n".format(*i)
            annotation_table.append(i)
        return [annotation_string, annotation_table]

    @staticmethod
    def add_cluster_table_to_database(cluster_table, db_path):
        '''get column names from Cluster object and create
        correspopnding table in database values from all
        clusters are filled to database'''
        column_name_and_type = []
        column_list = []

        # get all atribute names -> they are column names
        # in sqlite table, detect proper sqlite type
        def identity(x):
            return (x)

        for i in cluster_table[1]:
            t = type(cluster_table[1][i])
            if t == int:
                sqltype = "integer"
                convert = identity
            elif t == float:
                sqltype = "real"
                convert = identity
            elif t == bool:
                sqltype = "boolean"
                convert = bool
            else:
                sqltype = "text"
                convert = str
            column_name_and_type += ["[{}] {}".format(i, sqltype)]
            column_list += [tuple((i, convert))]
        header = ", ".join(column_name_and_type)
        db = sqlite3.connect(db_path)
        c = db.cursor()
        print("CREATE TABLE cluster_info ({})".format(header))
        c.execute("CREATE TABLE cluster_info ({})".format(header))
        # file data to cluster_table
        buffer = []
        for i in cluster_table:
            buffer.append(tuple('{}'.format(fun(i[j])) for j, fun in
                                column_list))
        wildcards = ",".join(["?"] * len(column_list))
        print(buffer)
        c.executemany("insert into cluster_info values ({})".format(wildcards),
                      buffer)
        db.commit()

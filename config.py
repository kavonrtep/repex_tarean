'''
All configuration for clustering
'''
import os
import tempfile
from math import exp
from collections import namedtuple
MAIN_DIR = os.path.dirname(os.path.realpath(__file__))
def add_base_path(base):
    '''automates generating absolute path in config'''
    def joined_path(p):
        '''create absolute path function '''
        return os.path.join(base, p)
    return joined_path

PATH = add_base_path(MAIN_DIR)

# clustering general settings
DIRECTORY_TREE = {'libdir': 'libdir',
                  'seqclust': 'seqclust',
                  'assembly': 'seqclust/small_clusters_assembly',
                  'blastx': 'seqclust/blastx',
                  'clustering': 'seqclust/clustering',
                  'clusters': 'seqclust/clustering/clusters',
                  'superclusters': 'seqclust/clustering/superclusters',
                  'mgblast': 'seqclust/mgblast',
                  'blastn': 'seqclust/blastn',
                  'prerun': 'seqclust/prerun',
                  'prerun_clusters': 'seqclust/prerun/clusters',
                  'sequences': 'seqclust/reads',
                  'custom_databases': 'seqclust/custom_databases'}

if "TEMP" in os.environ:
    DIRECTORY_TREE['TEMP'] = os.environ["TEMP"]
else:
    DIRECTORY_TREE['TEMP'] = tempfile.TemporaryDirectory().name

FILES = {'sample_db': DIRECTORY_TREE['TEMP'] + "/sample.db",
         'sample_fasta': DIRECTORY_TREE['prerun'] + "/sample.fasta",
         'prerun_cls_file' : DIRECTORY_TREE['prerun'] + "/sample_hitsort.cls",
         'filter_sequences_file' : DIRECTORY_TREE['prerun'] + "/filter_sequences.fasta",
         'sequences_db': DIRECTORY_TREE['TEMP'] + "/sequences.db",
         'sequences_fasta': DIRECTORY_TREE['sequences'] + "/reads.fasta",
         'hitsort': DIRECTORY_TREE['clustering'] + "/hitsort",
         'hitsort_db': DIRECTORY_TREE['TEMP'] + "/hitsort.db",
         'cls_file': DIRECTORY_TREE['clustering'] + "/hitsort.cls",
         'clusters_summary_csv': "CLUSTER_TABLE.csv",
         'profrep_classification_csv': "PROFREP_CLASSIFICATION_TEMPLATE.csv",
         'superclusters_csv_summary': "SUPERCLUSTER_TABLE.csv",
         'comparative_analysis_counts_csv': "COMPARATIVE_ANALYSIS_COUNTS.csv",
         'clusters_info': ".clusters_info.csv",
         'tarean_report_html': "tarean_report.html",
         'cluster_report_html' : "cluster_report.html",
         'supercluster_report_html' : 'supercluster_report.html',
         'repeat_annotation_summary_rds' : 'repeat_annotation_summary.rds',
         'summarized_annotation_html' :'summarized_annotation.html',
         'main_report_html' : 'index.html',
         'TR_consensus_fasta': "TAREAN_consensus_rank_{}.fasta",
         'summary_histogram' : 'summary_histogram.png',
         'comparative_summary_map': 'comparative_summary.png',
         "how_to_cite" : "HOW_TO_CITE.html",
         'logfile' : "logfile.txt",
         'contigs' : "contigs.fasta",
         'filter_omitted' : DIRECTORY_TREE['sequences'] + "/removed_filtering_positive_reads.fasta",
         'filter_kept' : DIRECTORY_TREE['sequences'] + "/kept_filtering_positive_reads.fasta"
}


# include in output-  [source, destination]
INCLUDE = [
    [PATH("HOW_TO_CITE.html"), FILES["how_to_cite"]]
]

# this is attribute of path - not a file name!
FILES_TO_DISCARD_AT_CLEANUP = [
    'prerun', 'mgblast', 'blastn', "blastx",
    'hitsort', "repeat_annotation_summary_rds"
]

# relative links for html files
HTML_LINKS = {
    "CLUSTER_TO_SUPERCLUSTER" : "../../superclusters/dir_SC%04d/index.html",
    "SUPERCLUSTER_TO_CLUSTER" : "../../clusters/dir_CL%04d/index.html",
    "CLUSTER_TO_CLUSTER" : "../dir_CL%04d/index.html",
    "SUPERCLUSTER_TO_SUPERCLUSTER" : "../dir_SC%04d/index.html",
    "CLUSTER_TO_CLUSTER_TABLE" : "../../../../cluster_report.html",
    "SEPERCLUSTER_TO_CLUSTER_TABLE" : "../../../../cluster_report.html",
    "ROOT_TO_CLUSTER" : "seqclust/clustering/clusters/dir_CL%04d/index.html",
    "ROOT_TO_SUPERCLUSTER" : "seqclust/clustering/superclusters/dir_SC%04d/index.html",
    "ROOT_TO_TAREAN" : "seqclust/clustering/clusters/dir_CL%04d/tarean/report.html",
    "CLUSTER_TO_KMER_REPORT" : "tarean/report.html",
    "INDEX_TO_TAREAN": "tarean_report.html",
    "INDEX_TO_CLUSTER_REPORT": "cluster_report.html",
    "INDEX_TO_SUPERCLUSTER_REPORT" : "supercluster_report.html",
    "INDEX_TO_SUMMARIZED_ANNOTATION" : "summarized_annotation.html"
}


EMAX = 42.6  # define how many graph edges can be processed in 1Kb RAM
# MINIMUM_NUMBER_OF_INPUT_SEQUENCES = 5000
# FOR TESTING:
MINIMUM_NUMBER_OF_INPUT_SEQUENCES = 1000
MINIMUM_NUMBER_OF_READS_IN_CLUSTER = 20  # smaller clusters are not analyzed
MINIMUM_NUMBER_OF_READS_FOR_MERGING = 20  # smaller clusters will not be merged
MINIMUM_NUMBER_OF_SHARED_PAIRS_FOR_MERGING = 20  # min size of W param


NUMBER_OF_SEQUENCES_FOR_PRERUN_WITH_FILTERING = 50000
NUMBER_OF_SEQUENCES_FOR_PRERUN_WITHOUT_FILTERING = 20000
NUMBER_OF_SEQUENCES_FOR_PRERUN = NUMBER_OF_SEQUENCES_FOR_PRERUN_WITHOUT_FILTERING
CHUNK_SIZE = 20000


CLUSTER_EMAX = 2E7  # this parameter higle affect memory usage!

CLUSTER_VMAX = 40000
SUPERCLUSTER_THRESHOLD = 0.1
# Number of processors to use - it will be set at runtime
PROC = None
RSERVE_PORT = 6311

#some settings related to repeats annotation
ORF_THRESHOLD = 1200
PBS_THRESHOLD = 2
# threshold for rDNA detection - percentage of similarity hits
RDNA_THRESHOLD = 20
# threshold for contamination detection
CONTAMINATION_THRESHOLD = 10
# tandem ranks codes:
# 1 : putative tandem repeats - high confidence
# 2 : putative tandem repeats - low confidence
# 3 : potential LTR element
# 4 : rDNA
TANDEM_RANKS = [1, 2, 3, 4]
SKIP_CAP3_ASSEMBLY_TANDEM_RANKS = [1]
FILTER_MIN_PROP_THRESHOLD = 0.03   # this is minimal proportion of graph edges!
FILTER_MIN_SIZE_THRESHOLD = 1000  # minimal size of the cluster to be consider for filtering
FILTER_PROPORTION_OF_KEPT = 0.1

R = 'lib'
# external scripts
RSOURCE_tarean = PATH('lib/tarean/tarean.R')
RSOURCE_reporting = PATH('lib/reporting.R')
RSOURCE_create_annotation = PATH('lib/create_annotation.R')
LTR_DETECTION = PATH("lib/detect_LTR_insertion_sites.pl")

#PATH to DATABASES:
DNA_DATABASE = PATH("databases/dna_database_masked.fasta")
TRNA_DATABASE = PATH("databases/tRNA_database.fasta")
SATELLITE_MODEL = PATH("databases/satellite_model.rds")
LASTAL_PARAMS = PATH("databases/lastal_params")
# for testing
PROTEIN_DATABASE = None
CLASSIFICATION_HIERARCHY = None
CUSTOM_DNA_DATABASE = None

# when modifying this section check if makefile has most recent target for protein database
PROTEIN_DATABASE_DEFAULT = "VIRIDIPLANTAE3.0"
PROTEIN_DATABASE_OPTIONS = {'VIRIDIPLANTAE3.0' :
                            (PATH("databases/protein_database_viridiplantae_v3.0.fasta"),   # change according if you use custom protein database
                             PATH(
                                 "databases/classification_tree_viridiplantae_v3.0_2.rds"
                                 "")),  # classification schem - data.tree object

                            'VIRIDIPLANTAE2.2' :
                            (PATH("databases/protein_database_viridiplantae_v2.2.fasta"),   # change according if you use custom protein database
                             PATH("databases/classification_viridiplantae_tree_2.rds")),  # classification schem - data.tree object
                            'METAZOA2.0' :
                            (PATH("databases/protein_database_metazoa_v2.fasta"),
                             # change according if you use custom protein database
                             PATH("databases/classification_tree_metazoa_v2_2.rds")),
                            # classification schem - data.tree object
                            'METAZOA3.0' :
                            (PATH("databases/protein_database_metazoa_v3.fasta"), # change according if you use custom protein database
                             PATH("databases/classification_tree_metazoa_v3_2.rds"))   # classification schem - data.tree object
}
# if you change PROTEIN_DATABASE_OPTIONS, do not forget to use 'makeblastdb' build blast database
# and 'diamond makedb' to build diamond database

# PATH to binaries
LOUVAIN = PATH("louvain")
BINARIES = PATH("bin")

CAP3_PATTERNS_REPLACE = {"{}.{}.contigs" : [">Contig", ">{}Contig"],
                         "{}.{}.aln" : ["* Contig", "* {}Contig"],
                         "{}.{}.contigs.qual" : [">Contig", ">{}Contig"],
                         "{}.{}.ace" : ["CO Contig", "CO {}Contig"],
                         "{}.{}.singlets" : None,
                         "{}.{}.info" : None,
                         "{}.{}.contigs.links" : None}

CAP3_FILES_MAPPING = {"{}.{}.contigs" : "small_clusters.fasta",
                      "{}.{}.contigs.qual" : "small_clusters.contigs.qual",
                      "{}.{}.aln" : "small_clusters.aln",
                      "{}.{}.ace" : "small_clusters.ace",
                      "{}.{}.singlets" : None,
                      "{}.{}.info" : None,
                      "{}.{}.contigs.links" : None}

CAP3_FILENAMES = list(CAP3_PATTERNS_REPLACE.keys())

CAP3_FILES_GOODNAMES = {"{}.{}.contigs" : "contigs.fasta",
                        "{}.{}.contigs.qual" : "contigs.qual",
                        "{}.{}.aln" : "contigs.aln",
                        "{}.{}.ace" : "contigs.ace",
                        "{}.{}.singlets" : "singlets.fasta",
                        "{}.{}.info" : "assembly.info",
                        "{}.{}.contigs.links" : "contigs.links"}




CAP3_PARAMS = " -p 80 -o 40 " ## not implemented yet

LTR_DETECTION_FILES = {'ADJ': 'LTR_info.ADJ',
                       'LTR': 'LTR_info.LTR',
                       'PBS_BLAST': 'LTR_info.with_PBS_blast.csv',
                       'BASE' : 'LTR_info'}

# options for all-2-all search and annotations
FilteringThreshod = namedtuple("FilteringThreshold",
                               "min_lcov min_pid min_ovl min_scov evalue")
AnnotationParams = namedtuple("AnnotationParams", "blastn blastx blastn_trna")
Option = namedtuple('Options',
                    ('name database all2all_search_params filtering_threshold '
                     'filter_self_hits legacy_database lastdb  annotation_search_params'))

# protein domain search options:
DIAMOND = {
    'args': ' -p {max_proc} --max-target-seqs 1 --min-score 30 --freq-sd 1000 --more-sensitive',
    'output_columns' : "qseqid sseqid qlen slen length ppos bitscore",
    'column_types' : [str, str, float, float, float, float, float],
    'program': 'diamond blastx',
    'filter_function' : lambda x: x.bitscore >= 30,
    'parallelize' : False
}
BLASTX_W3 = {
    'args': ' -num_alignments 1 -word_size 2 -evalue 0.01 ',
    'output_columns' : "qseqid sseqid qlen slen length ppos bitscore",
    'column_types' : [str, str, float, float, float, float, float],
    'program': 'blastx',
    'filter_function' : lambda x: x.bitscore >= 33
}
BLASTX_W2 = BLASTX_W3
BLASTX_W2['args'] = ' -num_alignments 1 -word_size 3 -evalue 0.01 '


ARGS = None

ILLUMINA = Option(
    name="illumina",
    database='blastdb_legacy',
    all2all_search_params=('mgblast -p 75 -W18 -UT -X40 -KT -JF  -F '
                           '"m D" -v100000000 -b100000000'
                           ' -D4 -C 30 -H 30 -i {query} -d {blastdb}'),
    filtering_threshold=FilteringThreshod(55, 90, 0, 0, 1),
    filter_self_hits=False,
    legacy_database=True,
    lastdb=False,
    annotation_search_params=AnnotationParams(
        blastn={
            'args': ' -task blastn  -num_alignments 1 -evalue 0.01 ',
            'output_columns' : "qseqid sseqid qlen slen length ppos bitscore",
            'column_types' : [str, str, float, float, float, float, float],
            'program': 'blastn',
            'filter_function' : lambda x: x.length > 30 and x.bitscore > 60

        },
        blastx=BLASTX_W3,
        blastn_trna={
            'args': ' -task blastn  -num_alignments 1 -word_size 7',
            'output_columns' : "qseqid sseqid qlen slen length ppos bitscore",
            'column_types' : [str, str, float, float, float, float, float],
            'program': 'blastn',
            'filter_function' : lambda x: x.length > 18 and x.bitscore > 60
        }
    )
)

ILLUMINA_DUST_OFF = ILLUMINA._replace(
    all2all_search_params=('mgblast -p 75 -W18 -UT -X40 -KT -JF  -F '
                           'F -v100000000 -b100000000'
                           ' -D4 -C 30 -H 30 -i {query} -d {blastdb}'),
)


ILLUMINA_SENSITIVE_MGBLAST = ILLUMINA._replace(
    all2all_search_params=('mgblast -p 75 -W8 -UT -X40 -KT -JF  -F '
                           'F -v100000000 -b100000000'
                           ' -D4 -C 30 -H 30 -i {query} -d {blastdb}'),
    filtering_threshold=FilteringThreshod(55, 80, 0, 0, 1),
)



ILLUMINA_SENSITIVE_BLASTPLUS = ILLUMINA._replace(
    name="illumina_short",
    database="blastdb",
    all2all_search_params=('blastplus_wrapper.py blastn -task blastn  -word_size 6 -perc_identity 80 -max_hsps 100000000 -max_target_seqs 100000000 -xdrop_gap 40 -xdrop_ungap 40'
                           ' -outfmt \\"6 qseqid qlen qstart qend sseqid slen sstart send pident score evalue sstrand\\" -query {query} -db {blastdb}'),
    filtering_threshold=FilteringThreshod(55, 80, 0, 0, 0.1)
)


OXFORD_NANOPORE = Option(
    name="oxford_nanopore",
    database='lastdb',
    all2all_search_params=('last_wrapper.py -f blasttab+ -P1 '
                           ' -m 700 -p {} '
                           ' {{blastdb}} {{query}}  ').format(LASTAL_PARAMS),
    filtering_threshold=FilteringThreshod(55, 50, 0, 0, 0.01),
    filter_self_hits=True,
    legacy_database=False,
    lastdb=True,
    annotation_search_params=AnnotationParams(
        blastn={
            'args': ' -task blastn  -num_alignments 1 -evalue 0.01 -word_size 11',
            'output_columns' : "qseqid sseqid qlen slen length ppos bitscore",
            'column_types' : [str, str, float, float, float, float, float],
            'program': 'blastn',
            'filter_function' : lambda x: x.length > 30 and x.bitscore > 50
        },
        blastx={
            'args': ' -num_alignments 1 -word_size 2 -evalue 0.1',
            'output_columns' : "qseqid sseqid qlen slen length ppos bitscore",
            'column_types' : [str, str, float, float, float, float, float],
            'program': 'blastx',
            'filter_function' : lambda x: x.bitscore >= 30
        },
        blastn_trna={
            'args': ' -task blastn  -num_alignments 1 -word_size 7',
            'output_columns' : "qseqid sseqid qlen slen length ppos bitscore",
            'column_types' : [str, str, float, float, float, float, float],
            'program': 'blastn',
            'filter_function' : lambda x: x.length > 18 and x.length > 60
        }
    )
)

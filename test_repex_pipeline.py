#!/usr/bin/env python3
'''
Basic Tarean and RepeatExplorer tests
'''
import subprocess
import tempfile
import unittest
import os
import shutil

def check_for_missing_files(directory, file_list):
    ''' check if files exists in the directory '''
    missing_files = []
    for f in file_list:
        path = os.path.join(directory, f)
        if os.path.exists(path):
            continue
        else:
            missing_files.append(f)
    return missing_files


class TestBasic(unittest.TestCase):
    ''' basic repex-tarean testcase '''
    EXECUTABLE = "./seqclust"

    # file lists to check
    FILE_LIST_BASIC = [
        "./seqclust/clustering/clusters/dir_CL0001/hitsort_part.csv",
        "./seqclust/clustering/clusters/dir_CL0001/reads.fasta",
        "./seqclust/clustering/clusters/dir_CL0001/reads_selection.fasta",
        "./seqclust/clustering/clusters/dir_CL0001/dna_database_annotation.csv",
        "./seqclust/clustering/clusters/dir_CL0001/graph_layout.GL",
        "./seqclust/clustering/clusters/dir_CL0001/graph_layout.png",
        "./seqclust/clustering/clusters/dir_CL0001/graph_layout_tmb.png",
        "./seqclust/clustering/clusters/dir_CL0001/graph_layout_directed.RData",
        "./logfile.txt", "./style1.css", "./documentation.html",
        "./tarean_report.html", "./cluster_report.html",
        "./summary_histogram.png", "./index.html", "./sequences.db",
        "./hitsort.db", "./TAREAN_consensus_rank_1.fasta",
        "./TAREAN_consensus_rank_2.fasta", "./TAREAN_consensus_rank_3.fasta",
        "./TAREAN_consensus_rank_4.fasta", "./seqclust/clustering/hitsort",
        "./seqclust/clustering/hitsort.cls"
    ]
    FILE_LIST_ASSEMBLY = [
        "./seqclust/small_clusters_assembly/small_clusters.aln",
        "./seqclust/small_clusters_assembly/small_clusters.ace",
        "./seqclust/small_clusters_assembly/small_clusters.fasta"
    ]
    FILE_LIST_FILTERING = ["./seqclust/prerun/filter_sequences.fasta"]
    FILE_LIST_COMPARATIVE = ["COMPARATIVE_ANALYSIS_COUNTS.csv"]
    FILE_LIST_CUSTOM_DATABASE = [
        "./seqclust/custom_databases/extra_database",
        "./seqclust/clustering/clusters/dir_CL0001/custom_db_extra_database_annotation.csv"
    ]
    def setUp(self):
        pass

    # helper function
    def tarean_run(self, cmd_options, file_list):
        ''' Basic taren run '''
        # output goes to tmp directory
        tmpdir = tempfile.mkdtemp()
        logfile = tempfile.NamedTemporaryFile(delete=False)
        print("\n------------------------------------------------------")
        print("Temp files:")
        print("   tmpdir : ", tmpdir)
        print("  logfile : ", logfile.name)
        print("------------------------------------------------------")
        print([self.EXECUTABLE] + ['-l', logfile.name, '-v', tmpdir] + cmd_options)
        p = subprocess.Popen(
            args=[self.EXECUTABLE] + ['-l', logfile.name, '-v', tmpdir
                                     ] + cmd_options)
        p.wait()
        status = p.returncode
        missing_files = check_for_missing_files(directory=tmpdir,
                                                file_list=file_list)
        if status:
            # print log file
            print("Non zero exit status!")
            with open(logfile.name) as f:
                print(f.read())

        self.assertEqual(status, 0)
        self.assertEqual(
            len(missing_files),
            0,
            msg="\n missing files: \n" + "\n".join(missing_files))
        shutil.rmtree(tmpdir)
        os.remove(logfile.name)


    def test_help(self):
        '''Test if help option works '''
        p = subprocess.Popen(args=[self.EXECUTABLE, "-h"],
                             stdout=subprocess.PIPE)
        output = str(p.stdout.readlines())
        p.stdout.close()
        p.wait()
        status = p.returncode
        self.assertRegex(output, "usage")
        self.assertRegex(output, "optional arguments")
        self.assertEqual(status, 0)

    def test_basic_no_merging_tarean(self):
        ''' Basic taren run '''
        cmd_options = ['-t', '-p', '-s', '6000', 'test_data/LAS_paired_10k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC)

    def test_basic_with_merging_tarean(self):
        ''' Basic taren run '''
        cmd_options = ['-t', '-p', '-M', '0.2', '-s', '6000',
                       'test_data/LAS_paired_10k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC)


    def test_basic_with_merging_tarean_dust_off(self):
        ''' Basic taren run '''
        cmd_options = ['-t', '-p', '-M', '0.2', '-s', '6000', "-opt", "ILLUMINA_DUST_OFF",
                       'test_data/LAS_paired_10k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC)

    def test_long_with_merging_tarean(self):
        '''Using more data with tarean'''
        cmd_options = ['-t', '-p', '-M', '0.1', '-m', '0.01',
                       'test_data/LAS_paired_25k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC)

    def test_long_with_merging2_tarean(self):
        '''Using more data with tarean 300k reads'''
        cmd_options = ['-t', '-p', '-M', '0.1', '-m', '0.01',
                       'test_data/LAS_paired_300k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC)

    def test_short_comparative_re(self):
        '''comparative analysis, two species, small run'''
        cmd_options = ['-P','3', '-p', '-m', '0.01',
                       'test_data/sequences_comparative.fasta']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC + self.FILE_LIST_COMPARATIVE)

    # REPEATEXPLORER - full runs
    def test_basic_no_merging_re(self):
        ''' Basic taren run '''
        cmd_options = ['-p', '-s', '6000', 'test_data/LAS_paired_10k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC + self.FILE_LIST_ASSEMBLY)

    def test_basic_no_merging_re_diamond(self):
        ''' Basic taren run '''
        cmd_options = ['-p', '-s', '6000','-D','DIAMOND', 'test_data/LAS_paired_10k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC + self.FILE_LIST_ASSEMBLY)



    def test_basic_with_merging_re(self):
        ''' Basic taren run '''
        cmd_options = ['-p', '-M', '0.2', '-s', '6000',
                       'test_data/LAS_paired_10k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC + self.FILE_LIST_ASSEMBLY)

    def test_long_with_merging_re(self):
        '''Using more data with tarean'''
        cmd_options = ['-p', '-M', '0.1', '-m', '0.01',
                       'test_data/LAS_paired_25k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC + self.FILE_LIST_ASSEMBLY)

    def test_long_with_merging_re_diamond(self):
        '''Using more data with tarean and using diamond'''
        cmd_options = ['-p', '-M', '0.1', '-m', '0.01','-D','DIAMOND',
                       'test_data/LAS_paired_25k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC + self.FILE_LIST_ASSEMBLY)

    def test_long_with_merging2_re(self):
        '''Using more data with tarean 300k reads'''
        cmd_options = ['-p', '-M', '0.1', '-m', '0.01',
                       'test_data/LAS_paired_300k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC + self.FILE_LIST_ASSEMBLY)

    def test_long_with_merging_and_filtering_re(self):
        '''Using more data with tarean, test of automatic filtering'''
        cmd_options = ['-A', '-p', '-M', '0.2', '-m', '0.01',
                       'test_data/ceu_200k.fasta']
        self.tarean_run(
            cmd_options,
            file_list=self.FILE_LIST_BASIC + self.FILE_LIST_FILTERING + self.FILE_LIST_ASSEMBLY)

    def test_custom_database_re(self):
        ''' Basic taren run '''
        cmd_options = ['-p', '-d', 'test_data/extra_database', 'extra_database', 'test_data/LAS_paired_10k.fas']
        self.tarean_run(cmd_options, file_list=self.FILE_LIST_BASIC + self.FILE_LIST_CUSTOM_DATABASE)

    def tearDown(self):
        pass


SHORT_TASK_NAME_LIST_TAREAN = ['test_help', 'test_basic_no_merging_tarean',
                               'test_basic_with_merging_tarean',
                               'test_basic_with_merging_tarean_dust_off']
LONG_TASK_NAME_LIST_TAREAN = ['test_long_with_merging_tarean',
                              'test_long_with_merging2_tarean']
SHORT_TASK_NAME_LIST_RE = ['test_basic_no_merging_re',
                           'test_basic_with_merging_re',
                           'test_basic_no_merging_re_diamond']
LONG_TASK_NAME_LIST_RE = ['test_long_with_merging_re',
                          'test_long_with_merging2_re',
                          'test_long_with_merging_and_filtering_re',
                          'test_long_with_merging_re_diamond']

COMPARATIVE_LIST = ['test_short_comparative_re']
CUSTOM_DATABASE_LIST = ['test_short_custom_database']

# Test suites:
SHORT_TAREAN_SUITE = unittest.TestSuite([TestBasic(i)
                                   for i in SHORT_TASK_NAME_LIST_TAREAN])
LONG_TAREAN_SUITE = unittest.TestSuite([TestBasic(i)
                                  for i in LONG_TASK_NAME_LIST_TAREAN])
COMPARATIVE_SUITE = unittest.TestSuite([TestBasic(i) for i in COMPARATIVE_LIST])
CUSTOM_DB_SUITE = unittest.TestSuite([TestBasic('test_custom_database_re')])

SHORT_RE_SUITE = unittest.TestSuite([TestBasic(i) for i in SHORT_TASK_NAME_LIST_RE])
LONG_RE_SUITE = unittest.TestSuite([TestBasic(i) for i in LONG_TASK_NAME_LIST_RE])

SHORT_SUITE = unittest.TestSuite([SHORT_RE_SUITE, SHORT_TAREAN_SUITE,
                                  COMPARATIVE_SUITE, CUSTOM_DB_SUITE])

LONG_LONG = unittest.TestSuite([LONG_RE_SUITE, LONG_TAREAN_SUITE])

# for single test tesing
if __name__ == '__main__':
    unittest.main(verbosity=2)

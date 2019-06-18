from IndexHoppingLib import *
from ErrorLib import *

from unittest import TestCase

import os
current_directory = os.getcwd()

class Test_parse_sample_sheet(TestCase):
    def setUp(self):
        self.test_sample_sheet_path = current_directory + "/test_files/TestSampleSheet.csv"
        self.expected_sample_sheet_info = {
            '17-ATZ02': {'Sample_Project': 'TEST-MSTR_TESTBATCH1_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH1'},
            '18-BNHL2': {'Sample_Project': 'TEST-MSTR_TESTBATCH2_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH2'},
            '18-BNHL9': {'Sample_Project': 'TEST-MSTR_TESTBATCH3_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH3'}}

    def test_parse_sample_sheet_happy(self):
        sample_sheet_info = parse_sample_sheet(self.test_sample_sheet_path)
        self.assertEqual(self.expected_sample_sheet_info, sample_sheet_info, "failed sample_sheet_info")


class Test_import_file(TestCase):
    def setUp(self):
        self.test_json_file = open(current_directory + "/test_files/Test_Stats.json", 'r')
        self.expected_flowcell_id = "TEST-FLOWCELL"
        self.expected_run_dir_base = "TEST_RUN_DIR_BASE"
        self.expected_conversion_results = [
            {'LaneNumber': 1, 'TotalClustersRaw': 610406, 'TotalClustersPF': 593007, 'Yield': 2965035, 'DemuxResults':
                [{'SampleId': '17-ATZ02', 'SampleName': '17-ATZ02', 'IndexMetrics':
                    [{'IndexSequence': 'AAAA+BBBB', 'MismatchCounts': {'0': 62916, '1': 1909}}],
                  'NumberReads': 1000, 'Yield': 324125,
                  'ReadMetrics': [{'ReadNumber': 1, 'Yield': 324125, 'YieldQ30': 305523, 'QualityScoreSum': 10740551,
                                   'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL2', 'SampleName': '18-BNHL2',
                  'IndexMetrics': [{'IndexSequence': 'CCCC+DDDD', 'MismatchCounts': {'0': 68117, '1': 1994}}],
                  'NumberReads': 1000, 'Yield': 350555,
                  'ReadMetrics': [{'ReadNumber': 1, 'Yield': 350555, 'YieldQ30': 327403, 'QualityScoreSum': 11717334,
                                   'TrimmedBases': 0}]}, {'SampleId': '18-BNHL9', 'SampleName': '18-BNHL9',
                                                          'IndexMetrics': [{'IndexSequence': 'EEEE+FFFF',
                                                                            'MismatchCounts': {'0': 1000, '1': 1865}}],
                                                          'NumberReads': 1000, 'Yield': 234965,
                                                          'ReadMetrics': [
                                                              {'ReadNumber': 1, 'Yield': 234965, 'YieldQ30': 207596,
                                                               'QualityScoreSum': 7638772, 'TrimmedBases': 0}]}],
             'Undetermined': {'NumberReads': 69538, 'Yield': 347690,
                              'ReadMetrics': [
                                  {'ReadNumber': 1, 'Yield': 347690, 'YieldQ30': 297761, 'QualityScoreSum': 11203807,
                                   'TrimmedBases': 0}]}},
            {'LaneNumber': 2, 'TotalClustersRaw': 610406, 'TotalClustersPF': 593007, 'Yield': 2965035, 'DemuxResults':
                [{'SampleId': '17-ATZ02', 'SampleName': '17-ATZ02', 'IndexMetrics':
                    [{'IndexSequence': 'AAAA+BBBB', 'MismatchCounts': {'0': 62916, '1': 1909}}],
                  'NumberReads': 1000, 'Yield': 324125, 'ReadMetrics':
                      [{'ReadNumber': 1, 'Yield': 324125, 'YieldQ30': 305523, 'QualityScoreSum': 10740551,
                        'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL2', 'SampleName': '18-BNHL2',
                  'IndexMetrics': [{'IndexSequence': 'CCCC+DDDD', 'MismatchCounts':
                      {'0': 68117, '1': 1994}}], 'NumberReads': 1000, 'Yield': 350555, 'ReadMetrics':
                      [{'ReadNumber': 1, 'Yield': 350555, 'YieldQ30': 327403, 'QualityScoreSum': 11717334,
                        'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL9', 'SampleName': '18-BNHL9',
                  'IndexMetrics': [{'IndexSequence': 'EEEE+FFFF', 'MismatchCounts': {'0': 1000, '1': 1865}}],
                  'NumberReads': 1000, 'Yield': 234965, 'ReadMetrics': [
                     {'ReadNumber': 1, 'Yield': 234965, 'YieldQ30': 207596, 'QualityScoreSum': 7638772,
                      'TrimmedBases': 0}]}],
             'Undetermined': {'NumberReads': 69538, 'Yield': 347690, 'ReadMetrics': [
                 {'ReadNumber': 1, 'Yield': 347690, 'YieldQ30': 297761, 'QualityScoreSum': 11203807,
                  'TrimmedBases': 0}]}}]
        self.expected_num_of_lanes = 2
        self.expected_num_of_samples = 3
        self.expected_unknown_barcodes = [{'Lane': 1, 'Barcodes':
            {'AAAB+DDDD': 200, 'AAAB+BBBB': 300, 'AAAA+DDDD': 100, 'CCCC+BBBB': 100, 'CCCB+FFFZ': 30, 'AAAA+ZZZZ': 250,
             'ZZZZ+DDDD': 100, 'ABBA+BBBB': 50, 'ZZZZ+ZZZZ': 500, 'GGGG+FFDF': 40, 'BBBB+AAAA': 1000}},
                                          {'Lane': 2, 'Barcodes':
                                              {'AAAB+DDDD': 1, 'AAAB+BBBB': 1, 'AAAA+DDDD': 1, 'CCCC+BBBB': 1,
                                               'CCCB+FFFZ': 1, 'AAAA+ZZZZ': 1, 'ZZZZ+DDDD': 1, 'ABBA+BBBB': 1,
                                               'ZZZZ+ZZZZ': 1, 'GGGG+FFDF': 1, 'ABCD+DDDE': 1, 'BBBB+AAAA': 1}}]
        self.expected_sample_list = ['17-ATZ02', '18-BNHL2', '18-BNHL9']

    def test_import_file_happy(self):
        flowcell_id, run_dir_base, conversion_results, num_of_lanes, num_of_samples, unknown_barcodes, sample_list = import_file(
            self.test_json_file)
        self.test_json_file.close()
        self.assertEqual(self.expected_flowcell_id, flowcell_id, "failed flowcell_id")
        self.assertEqual(self.expected_run_dir_base, run_dir_base, "failed run_dir_base")
        self.assertEqual(self.expected_conversion_results, conversion_results, "failed conversion_results")
        self.assertEqual(self.expected_num_of_lanes, num_of_lanes, "failed num_of_lanes")
        self.assertEqual(self.expected_num_of_samples, num_of_samples, "failed num_of_samples")
        self.assertEqual(self.expected_unknown_barcodes, unknown_barcodes, "failed unknown_barcodes")
        self.assertEqual(self.expected_sample_list, sample_list, "failed sample_list")


class Test_capture_index_sequence(TestCase):
    def setUp(self):
        self.conversion_results = [
            {'LaneNumber': 1, 'TotalClustersRaw': 610406, 'TotalClustersPF': 593007, 'Yield': 2965035, 'DemuxResults':
                [{'SampleId': '17-ATZ02', 'SampleName': '17-ATZ02', 'IndexMetrics':
                    [{'IndexSequence': 'AAAA+BBBB', 'MismatchCounts': {'0': 62916, '1': 1909}}],
                  'NumberReads': 1000, 'Yield': 324125,
                  'ReadMetrics': [{'ReadNumber': 1, 'Yield': 324125, 'YieldQ30': 305523, 'QualityScoreSum': 10740551,
                                   'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL2', 'SampleName': '18-BNHL2',
                  'IndexMetrics': [{'IndexSequence': 'CCCC+DDDD', 'MismatchCounts': {'0': 68117, '1': 1994}}],
                  'NumberReads': 1000, 'Yield': 350555,
                  'ReadMetrics': [{'ReadNumber': 1, 'Yield': 350555, 'YieldQ30': 327403, 'QualityScoreSum': 11717334,
                                   'TrimmedBases': 0}]}, {'SampleId': '18-BNHL9', 'SampleName': '18-BNHL9',
                                                          'IndexMetrics': [{'IndexSequence': 'EEEE+FFFF',
                                                                            'MismatchCounts': {'0': 1000, '1': 1865}}],
                                                          'NumberReads': 1000, 'Yield': 234965,
                                                          'ReadMetrics': [
                                                              {'ReadNumber': 1, 'Yield': 234965, 'YieldQ30': 207596,
                                                               'QualityScoreSum': 7638772, 'TrimmedBases': 0}]}],
             'Undetermined': {'NumberReads': 69538, 'Yield': 347690,
                              'ReadMetrics': [
                                  {'ReadNumber': 1, 'Yield': 347690, 'YieldQ30': 297761, 'QualityScoreSum': 11203807,
                                   'TrimmedBases': 0}]}},
            {'LaneNumber': 2, 'TotalClustersRaw': 610406, 'TotalClustersPF': 593007, 'Yield': 2965035, 'DemuxResults':
                [{'SampleId': '17-ATZ02', 'SampleName': '17-ATZ02', 'IndexMetrics':
                    [{'IndexSequence': 'AAAA+BBBB', 'MismatchCounts': {'0': 62916, '1': 1909}}],
                  'NumberReads': 1000, 'Yield': 324125, 'ReadMetrics':
                      [{'ReadNumber': 1, 'Yield': 324125, 'YieldQ30': 305523, 'QualityScoreSum': 10740551,
                        'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL2', 'SampleName': '18-BNHL2',
                  'IndexMetrics': [{'IndexSequence': 'CCCC+DDDD', 'MismatchCounts':
                      {'0': 68117, '1': 1994}}], 'NumberReads': 1000, 'Yield': 350555, 'ReadMetrics':
                      [{'ReadNumber': 1, 'Yield': 350555, 'YieldQ30': 327403, 'QualityScoreSum': 11717334,
                        'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL9', 'SampleName': '18-BNHL9',
                  'IndexMetrics': [{'IndexSequence': 'EEEE+FFFF', 'MismatchCounts': {'0': 1000, '1': 1865}}],
                  'NumberReads': 1000, 'Yield': 234965, 'ReadMetrics': [
                     {'ReadNumber': 1, 'Yield': 234965, 'YieldQ30': 207596, 'QualityScoreSum': 7638772,
                      'TrimmedBases': 0}]}],
             'Undetermined': {'NumberReads': 69538, 'Yield': 347690, 'ReadMetrics': [
                 {'ReadNumber': 1, 'Yield': 347690, 'YieldQ30': 297761, 'QualityScoreSum': 11203807,
                  'TrimmedBases': 0}]}}]
        self.conversion_results_no_dual_index = [
            {'LaneNumber': 1, 'TotalClustersRaw': 610406, 'TotalClustersPF': 593007, 'Yield': 2965035, 'DemuxResults':
                [{'SampleId': '17-ATZ02', 'SampleName': '17-ATZ02', 'IndexMetrics':
                    [{'IndexSequence': 'AAAA+BBBB', 'MismatchCounts': {'0': 62916, '1': 1909}}],
                  'NumberReads': 1000, 'Yield': 324125,
                  'ReadMetrics': [{'ReadNumber': 1, 'Yield': 324125, 'YieldQ30': 305523, 'QualityScoreSum': 10740551,
                                   'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL2', 'SampleName': '18-BNHL2',
                  'IndexMetrics': [{'IndexSequence': 'CCCC', 'MismatchCounts': {'0': 68117, '1': 1994}}],
                  'NumberReads': 1000, 'Yield': 350555,
                  'ReadMetrics': [{'ReadNumber': 1, 'Yield': 350555, 'YieldQ30': 327403, 'QualityScoreSum': 11717334,
                                   'TrimmedBases': 0}]}, {'SampleId': '18-BNHL9', 'SampleName': '18-BNHL9',
                                                          'IndexMetrics': [{'IndexSequence': 'EEEE+FFFF',
                                                                            'MismatchCounts': {'0': 1000, '1': 1865}}],
                                                          'NumberReads': 1000, 'Yield': 234965,
                                                          'ReadMetrics': [
                                                              {'ReadNumber': 1, 'Yield': 234965, 'YieldQ30': 207596,
                                                               'QualityScoreSum': 7638772, 'TrimmedBases': 0}]}],
             'Undetermined': {'NumberReads': 69538, 'Yield': 347690,
                              'ReadMetrics': [
                                  {'ReadNumber': 1, 'Yield': 347690, 'YieldQ30': 297761, 'QualityScoreSum': 11203807,
                                   'TrimmedBases': 0}]}},
            {'LaneNumber': 2, 'TotalClustersRaw': 610406, 'TotalClustersPF': 593007, 'Yield': 2965035, 'DemuxResults':
                [{'SampleId': '17-ATZ02', 'SampleName': '17-ATZ02', 'IndexMetrics':
                    [{'IndexSequence': 'AAAA+BBBB', 'MismatchCounts': {'0': 62916, '1': 1909}}],
                  'NumberReads': 1000, 'Yield': 324125, 'ReadMetrics':
                      [{'ReadNumber': 1, 'Yield': 324125, 'YieldQ30': 305523, 'QualityScoreSum': 10740551,
                        'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL2', 'SampleName': '18-BNHL2',
                  'IndexMetrics': [{'IndexSequence': 'CCCC', 'MismatchCounts':
                      {'0': 68117, '1': 1994}}], 'NumberReads': 1000, 'Yield': 350555, 'ReadMetrics':
                      [{'ReadNumber': 1, 'Yield': 350555, 'YieldQ30': 327403, 'QualityScoreSum': 11717334,
                        'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL9', 'SampleName': '18-BNHL9',
                  'IndexMetrics': [{'IndexSequence': 'EEEE+FFFF', 'MismatchCounts': {'0': 1000, '1': 1865}}],
                  'NumberReads': 1000, 'Yield': 234965, 'ReadMetrics': [
                     {'ReadNumber': 1, 'Yield': 234965, 'YieldQ30': 207596, 'QualityScoreSum': 7638772,
                      'TrimmedBases': 0}]}],
             'Undetermined': {'NumberReads': 69538, 'Yield': 347690, 'ReadMetrics': [
                 {'ReadNumber': 1, 'Yield': 347690, 'YieldQ30': 297761, 'QualityScoreSum': 11203807,
                  'TrimmedBases': 0}]}}]
        self.expected_index_sequence = ['AAAA+BBBB', 'CCCC+DDDD', 'EEEE+FFFF']

    def test_capture_index_sequence_happy(self):
        index_sequence = capture_index_sequence(self.conversion_results)
        self.assertEqual(self.expected_index_sequence, index_sequence, "failed index_sequence")
    def test_capture_index_sequence_no_dual_index(self):
        with self.assertRaises(NoDualIndexException):
            capture_index_sequence(self.conversion_results_no_dual_index)


class Test_make_all_possible_index_combinations(TestCase):
    def setUp(self):
        self.index_sequence = ['AAAA+BBBB', 'CCCC+DDDD', 'EEEE+FFFF']
        self.index_sequence_not_unique1 = ['AAAA+BBBB', 'AAAA+DDDD', 'EEEE+FFFF']
        self.index_sequence_not_unique2 = ['AAAA+BBBB', 'CCCC+DDDD', 'EEEE+BBBB']
        self.expected_index1_sequence = ['AAAA', 'CCCC', 'EEEE']
        self.expected_index2_sequence = ['BBBB', 'DDDD', 'FFFF']
        self.expected_mismatch_index_sequence = ['AAAA+DDDD', 'AAAA+FFFF', 'CCCC+BBBB', 'CCCC+FFFF', 'EEEE+BBBB',
                                                 'EEEE+DDDD']

    def test_make_all_possible_index_combinations_happy(self):
        index1_sequence, index2_sequence, mismatch_index_sequence = make_all_possible_index_combinations(
            self.index_sequence)
        self.assertEqual(self.expected_index1_sequence, index1_sequence, "failed make_all_possible_index_combination: "
                                                                         "index1_sequence")
        self.assertEqual(self.expected_index2_sequence, index2_sequence, "failed make_all_possible_index_combination: "
                                                                         "index2_sequence")
        self.assertEqual(self.expected_mismatch_index_sequence, mismatch_index_sequence, "failed make_all_possible_"
                                                                                         "index_combination: mismatch_index_sequence")

    def test_make_all_possible_index_combinations_not_unique1(self):
        with self.assertRaises(NotUniqueDualIndexException):
            make_all_possible_index_combinations(self.index_sequence_not_unique1)

    def test_make_all_possible_index_combinations_not_unique2(self):
        with self.assertRaises(NotUniqueDualIndexException):
            make_all_possible_index_combinations(self.index_sequence_not_unique2)

class Test_compare_sequences(TestCase):
    def setUp(self):
        self.sequence1 = 'AAAA'
        self.sequence2 = 'AAAB'
        self.expected_similarity_score = 1

    def test_compare_sequences_happy(self):
        similarity_score = compare_sequences(self.sequence1, self.sequence2)
        self.assertEqual(self.expected_similarity_score, similarity_score, "failed compare_sequences: similarity score")


class Test_calculate_total_number_of_reads(TestCase):
    def setUp(self):
        self.conversion_results = [
            {'LaneNumber': 1, 'TotalClustersRaw': 610406, 'TotalClustersPF': 593007, 'Yield': 2965035, 'DemuxResults':
                [{'SampleId': '17-ATZ02', 'SampleName': '17-ATZ02', 'IndexMetrics':
                    [{'IndexSequence': 'AAAA+BBBB', 'MismatchCounts': {'0': 62916, '1': 1909}}],
                  'NumberReads': 1000, 'Yield': 324125,
                  'ReadMetrics': [{'ReadNumber': 1, 'Yield': 324125, 'YieldQ30': 305523, 'QualityScoreSum': 10740551,
                                   'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL2', 'SampleName': '18-BNHL2',
                  'IndexMetrics': [{'IndexSequence': 'CCCC+DDDD', 'MismatchCounts': {'0': 68117, '1': 1994}}],
                  'NumberReads': 1000, 'Yield': 350555,
                  'ReadMetrics': [{'ReadNumber': 1, 'Yield': 350555, 'YieldQ30': 327403, 'QualityScoreSum': 11717334,
                                   'TrimmedBases': 0}]}, {'SampleId': '18-BNHL9', 'SampleName': '18-BNHL9',
                                                          'IndexMetrics': [{'IndexSequence': 'EEEE+FFFF',
                                                                            'MismatchCounts': {'0': 1000, '1': 1865}}],
                                                          'NumberReads': 1000, 'Yield': 234965,
                                                          'ReadMetrics': [
                                                              {'ReadNumber': 1, 'Yield': 234965, 'YieldQ30': 207596,
                                                               'QualityScoreSum': 7638772, 'TrimmedBases': 0}]}],
             'Undetermined': {'NumberReads': 69538, 'Yield': 347690,
                              'ReadMetrics': [
                                  {'ReadNumber': 1, 'Yield': 347690, 'YieldQ30': 297761, 'QualityScoreSum': 11203807,
                                   'TrimmedBases': 0}]}},
            {'LaneNumber': 2, 'TotalClustersRaw': 610406, 'TotalClustersPF': 593007, 'Yield': 2965035, 'DemuxResults':
                [{'SampleId': '17-ATZ02', 'SampleName': '17-ATZ02', 'IndexMetrics':
                    [{'IndexSequence': 'AAAA+BBBB', 'MismatchCounts': {'0': 62916, '1': 1909}}],
                  'NumberReads': 1000, 'Yield': 324125, 'ReadMetrics':
                      [{'ReadNumber': 1, 'Yield': 324125, 'YieldQ30': 305523, 'QualityScoreSum': 10740551,
                        'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL2', 'SampleName': '18-BNHL2',
                  'IndexMetrics': [{'IndexSequence': 'CCCC+DDDD', 'MismatchCounts':
                      {'0': 68117, '1': 1994}}], 'NumberReads': 1000, 'Yield': 350555, 'ReadMetrics':
                      [{'ReadNumber': 1, 'Yield': 350555, 'YieldQ30': 327403, 'QualityScoreSum': 11717334,
                        'TrimmedBases': 0}]},
                 {'SampleId': '18-BNHL9', 'SampleName': '18-BNHL9',
                  'IndexMetrics': [{'IndexSequence': 'EEEE+FFFF', 'MismatchCounts': {'0': 1000, '1': 1865}}],
                  'NumberReads': 1000, 'Yield': 234965, 'ReadMetrics': [
                     {'ReadNumber': 1, 'Yield': 234965, 'YieldQ30': 207596, 'QualityScoreSum': 7638772,
                      'TrimmedBases': 0}]}],
             'Undetermined': {'NumberReads': 69538, 'Yield': 347690, 'ReadMetrics': [
                 {'ReadNumber': 1, 'Yield': 347690, 'YieldQ30': 297761, 'QualityScoreSum': 11203807,
                  'TrimmedBases': 0}]}}]
        self.num_of_lanes = 2
        self.num_of_samples = 3
        self.expected_total_number_of_reads = 6000

    def test_calculate_total_number_of_reads_happy(self):
        total_number_of_reads = calculate_total_number_of_reads(self.conversion_results, self.num_of_lanes,
                                                                self.num_of_samples)
        self.assertEqual(self.expected_total_number_of_reads, total_number_of_reads,
                         "failed calculate_total_number_of_reads: total_number_of_reads")


class Test_mismatched_reads(TestCase):
    def setUp(self):
        self.index1_sequence = ['AAAA', 'CCCC', 'EEEE']
        self.index2_sequence = ['BBBB', 'DDDD', 'FFFF']
        self.mismatch_index_sequence = ['AAAA+DDDD', 'AAAA+FFFF', 'CCCC+BBBB', 'CCCC+FFFF', 'EEEE+BBBB', 'EEEE+DDDD']
        self.num_of_lanes = 2
        self.unknown_barcodes = [{'Lane': 1, 'Barcodes': {'AAAB+DDDD': 200, 'AAAB+BBBB': 300, 'AAAA+DDDD': 100,
                                                          'CCCC+BBBB': 100, 'CCCB+FFFZ': 30, 'AAAA+ZZZZ': 250,
                                                          'ZZZZ+DDDD': 100,
                                                          'ABBA+BBBB': 50, 'ZZZZ+ZZZZ': 500, 'GGGG+FFDF': 40,
                                                          'BBBB+AAAA': 1000}},
                                 {'Lane': 2,
                                  'Barcodes': {'AAAB+DDDD': 1, 'AAAB+BBBB': 1, 'AAAA+DDDD': 1, 'CCCC+BBBB': 1,
                                               'CCCB+FFFZ': 1, 'AAAA+ZZZZ': 1, 'ZZZZ+DDDD': 1, 'ABBA+BBBB': 1,
                                               'ZZZZ+ZZZZ': 1, 'GGGG+FFDF': 1, 'ABCD+DDDE': 1, 'BBBB+AAAA': 1}}]
        self.expected_mismatch_index_dict = {'AAAA+DDDD': 101, 'CCCC+BBBB': 101}
        self.expected_similar_mismatch_index_dict = {'AAAB+DDDD': 201, 'CCCB+FFFZ': 31}
        self.expected_not_similar_mismatch_index_dict = {'AAAA+ZZZZ': 251, 'ZZZZ+DDDD': 101, 'ABBA+BBBB': 51,
                                                         'GGGG+FFDF': 41, 'ABCD+DDDE': 1}

    def test_mismatched_reads_happy(self):
        mismatch_index_dict, similar_mismatch_index_dict, not_similar_mismatch_index_dict = mismatched_reads(
            self.index1_sequence,
            self.index2_sequence,
            self.mismatch_index_sequence,
            self.num_of_lanes,
            self.unknown_barcodes)
        self.assertEqual(self.expected_mismatch_index_dict, mismatch_index_dict,
                         "failed mismatched_reads: mismatch_index_dict")
        self.assertEqual(self.expected_not_similar_mismatch_index_dict, not_similar_mismatch_index_dict,
                         "failed mismatched_reads: not_similar_mismatch_index_dict")
        self.assertEqual(self.expected_similar_mismatch_index_dict, similar_mismatch_index_dict,
                         "failed mismatched reads: similar_mismatch_index_dict")


class Test_index_hopping_percent(TestCase):
    def setUp(self):
        self.mismatch_index_dict = {'AAAA+DDDD': 101, 'CCCC+BBBB': 101}
        self.similar_mismatch_index_dict = {'AAAB+DDDD': 201, 'CCCB+FFFZ': 31}
        self.total_number_of_reads = 6000
        self.expected_index_hop_perecent = 7.23
        self.expected_total_number_of_mismatched_reads = 434

    def test_index_hopping_percent_happy(self):
        index_hop_percent, total_number_of_mismatched_reads = index_hopping_percent(self.mismatch_index_dict,
                                                                                    self.similar_mismatch_index_dict,
                                                                                    self.total_number_of_reads)
        self.assertAlmostEqual(self.expected_index_hop_perecent, index_hop_percent, places=2,
                               msg="failed index_hopping_percent: index_hop_percent")
        self.assertEqual(self.expected_total_number_of_mismatched_reads, total_number_of_mismatched_reads,
                         msg="failed index_hopping_percent: total_number_of_mismatched_reads")


class Test_index_jump_count(TestCase):
    def setUp(self):
        self.index1_sequence = ['AAAA', 'CCCC', 'EEEE']
        self.index2_sequence = ['BBBB', 'DDDD', 'FFFF']
        self.index_sequence = ['AAAA+BBBB', 'CCCC+DDDD', 'EEEE+FFFF']
        self.mismatch_index_dict = {'AAAA+DDDD': 101, 'CCCC+BBBB': 101}
        self.similar_mismatch_index_dict = {'AAAB+DDDD': 201, 'CCCB+FFFZ': 31}
        self.expected_index_jump_dict = {'AAAA+BBBB': 403, 'CCCC+DDDD': 434, 'EEEE+FFFF': 31}

    def test_index_jump_count_happy(self):
        index_jump_dict = index_jump_count(self.index1_sequence, self.index2_sequence, self.index_sequence,
                                           self.mismatch_index_dict, self.similar_mismatch_index_dict)
        self.assertEqual(self.expected_index_jump_dict, index_jump_dict, msg="index_jump_count_happy")


class Test_not_similar_jump_count(TestCase):
    def setUp(self):
        self.index1_sequence = ['AAAA', 'CCCC', 'EEEE']
        self.index2_sequence = ['BBBB', 'DDDD', 'FFFF']
        self.index_sequence = ['AAAA+BBBB', 'CCCC+DDDD', 'EEEE+FFFF']
        self.not_similar_mismatch_index_dict = {'AAAA+ZZZZ': 251, 'ZZZZ+DDDD': 101, 'ABBA+BBBB': 51,
                                                'GGGG+FFDF': 41, 'ABCD+DDDE': 1}
        self.expected_not_similar_index_association = {'AAAA+BBBB': {'AAAA+ZZZZ': 251, 'ABBA+BBBB': 51},
                                                       'CCCC+DDDD': {'ZZZZ+DDDD': 101, 'ABCD+DDDE': 1},
                                                       'EEEE+FFFF': {'GGGG+FFDF': 41}}
        self.expected_not_similar_jump_count_dict = {'AAAA+BBBB': 302, 'CCCC+DDDD': 102, 'EEEE+FFFF': 41}

    def test_not_similar_jump_count_happy(self):
        not_similar_index_association, not_similar_jump_count_dict = not_similar_jump_count(self.index1_sequence,
                                                                                            self.index2_sequence,
                                                                                            self.index_sequence,
                                                                                            self.not_similar_mismatch_index_dict)
        self.assertEqual(self.expected_not_similar_index_association, not_similar_index_association,
                         msg="not_similar_index_association")
        self.assertEqual(self.expected_not_similar_jump_count_dict, not_similar_jump_count_dict,
                         msg="not_similar_jump_count_dict")


class Test_validate_samples(TestCase):
    def setUp(self):
        self.sample_list = ['17-ATZ02', '18-BNHL2', '18-BNHL9']
        self.sample_sheet_info = {
            '17-ATZ02': {'Sample_Project': 'TEST-MSTR_TESTBATCH1_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH1'},
            '18-BNHL2': {'Sample_Project': 'TEST-MSTR_TESTBATCH2_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH2'},
            '18-BNHL9': {'Sample_Project': 'TEST-MSTR_TESTBATCH3_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH3'}}
        self.sample_sheet_info_wrong_sample = {
            'XXXX': {'Sample_Project': 'TEST-MSTR_TESTBATCH1_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                     'Batch_ID': 'TESTBATCH1'},
            '18-BNHL2': {'Sample_Project': 'TEST-MSTR_TESTBATCH2_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH2'},
            '18-BNHL9': {'Sample_Project': 'TEST-MSTR_TESTBATCH3_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH3'}}
        self.sample_sheet_info_wrong_number = {
            '17-ATZ02': {'Sample_Project': 'TEST-MSTR_TESTBATCH1_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                     'Batch_ID': 'TESTBATCH1'},
            '18-BNHL2': {'Sample_Project': 'TEST-MSTR_TESTBATCH2_TESTCHEMID_TESTSEQID_TESTFLOWCELLID',
                         'Batch_ID': 'TESTBATCH2'}}

    def test_validate_samples_error(self):
        with self.assertRaises(SampleNotFoundException):
            validate_samples(self.sample_list, self.sample_sheet_info_wrong_sample)

    def test_validate_samples_unequal_numbers(self):
        with self.assertRaises(SampleNumberMismatchException):
            validate_samples(self.sample_list, self.sample_sheet_info_wrong_number)

    def test_validate_samples_happy(self):
        validate_samples(self.sample_list, self.sample_sheet_info)

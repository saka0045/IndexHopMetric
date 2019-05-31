import ParseJson
from unittest import TestCase

from ParseJson import index_jump_count


class TestIndex_jump_count(TestCase):

    def setUp(self):
        self.index1_sequence = ['CAGCAGGTCA', 'ATAAGCGGAG', 'TAACTTGGTC', 'CGACCTTAAT', 'TCTAGCTGCG', 'AGATTAGCGT',
                           'AGGTAACAAC', 'ATGTTGCCAC']
        self.index_sequence = ['CAGCAGGTCA+ACCTAAGTGA', 'ATAAGCGGAG+GTAGAACAGA', 'TAACTTGGTC+ATTCACGACA',
                          'CGACCTTAAT+CCTAAGAGCA', 'TCTAGCTGCG+CTTACCTATA', 'AGATTAGCGT+TATAAGGCGA',
                          'AGGTAACAAC+GCATCGTGCA', 'ATGTTGCCAC+TGGCTTCTAA']
        self.index2_sequence = ['ACCTAAGTGA', 'GTAGAACAGA', 'ATTCACGACA', 'CCTAAGAGCA', 'CTTACCTATA', 'TATAAGGCGA',
                           'GCATCGTGCA', 'TGGCTTCTAA']
        self.mismatch_index_dict = {'AGATTAGCGT+CTTACCTATA': 300, 'CGACCTTAAT+TGGCTTCTAA': 80, 'CAGCAGGTCA+CTTACCTATA': 80,
                               'TAACTTGGTC+GTAGAACAGA': 60, 'TAACTTGGTC+CTTACCTATA': 60, 'AGGTAACAAC+TATAAGGCGA': 60,
                               'AGGTAACAAC+CTTACCTATA': 60, 'TAACTTGGTC+TGGCTTCTAA': 40, 'CGACCTTAAT+TATAAGGCGA': 40,
                               'CGACCTTAAT+CTTACCTATA': 40, 'ATAAGCGGAG+CCTAAGAGCA': 40, 'AGATTAGCGT+TGGCTTCTAA': 40,
                               'TCTAGCTGCG+TATAAGGCGA': 20, 'TCTAGCTGCG+GCATCGTGCA': 20, 'TCTAGCTGCG+CCTAAGAGCA': 20,
                               'TAACTTGGTC+ACCTAAGTGA': 20, 'CGACCTTAAT+GTAGAACAGA': 20, 'CGACCTTAAT+GCATCGTGCA': 20,
                               'CGACCTTAAT+ATTCACGACA': 20}
        self.similar_mismatch_index_dict = {'TTACTTGGTC+CTTACCTATA': 20, 'TGATTAGCGT+CTTACCTATA': 20,
                                       'TGACTTGGTC+TATAAGGCCA': 20, 'TGACCTTAAT+TATAAGGCGA': 20,
                                       'TCTAGCTGCT+TATAAGGCGA': 20, 'TCTAGCTGCG+TATACGGCGA': 20,
                                       'TCAAGCTGCG+TATAAGGCGA': 20, 'TATCTTGGTC+CTTCCCTATA': 20,
                                       'TAACTTGGTC+TATAACGCGA': 20, 'TAACTTGGTC+CTTCCCTATA': 20,
                                       'TAACTTGGTC+CTTACCTACA': 20, 'TAACTTGGTC+CTTACCGATA': 20,
                                       'CTGTTGCCAC+CTTACCTATA': 20, 'CGTCCTTAAT+CTTACCTATA': 20,
                                       'CGGCCTTAAT+GCATCGTGCA': 20, 'CGACCTTCAT+TGGCTTCTAA': 20,
                                       'CGACCTTAGT+TGGCTTCTAA': 20, 'CGACCTTAAT+GCATCGGGCA': 20,
                                       'CGACCTTAAT+GCATCGAGCA': 20, 'CGACCTTAAT+GCAACGTGCA': 20,
                                       'CGACCTTAAT+CTTCACGACA': 20, 'CGACCTTAAT+CTTACCTAAA': 20,
                                       'CGACCTTAAT+AACTAAGTGA': 20}
        self.expected_index_jump_dict = {'CAGCAGGTCA+ACCTAAGTGA': 120, 'ATAAGCGGAG+GTAGAACAGA': 120,
                                    'TAACTTGGTC+ATTCACGACA': 360, 'CGACCTTAAT+CCTAAGAGCA': 500,
                                    'TCTAGCTGCG+CTTACCTATA': 840, 'AGATTAGCGT+TATAAGGCGA': 600,
                                    'AGGTAACAAC+GCATCGTGCA': 240, 'ATGTTGCCAC+TGGCTTCTAA': 220}

    def test_index_jump_count_happy(self):
        index_jump_dict = index_jump_count(self.index1_sequence, self.index2_sequence, self.index_sequence, self.mismatch_index_dict, self.similar_mismatch_index_dict)
        self.assertEqual(self.expected_index_jump_dict, index_jump_dict, "Expected index dictionary")

    def test_index_jump_count_no_index1(self):
        with self.assertRaises(ParseJson.MissingSequence1Exception):
            index_jump_count([], self.index2_sequence, self.index_sequence, self.mismatch_index_dict, self.similar_mismatch_index_dict)



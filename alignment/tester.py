import dynamic_practice
import dag
import aligner
import unittest

class Tester(unittest.TestCase):
    known_coins = ((40, [1, 5, 10, 20, 25, 50], 2),)

    known_length = (([[1, 0, 2, 4, 3], [4, 6, 5, 2, 1], [4, 4, 5, 2, 1],
                      [5, 6, 8, 5, 3]],
                     [[3, 2, 4, 0], [3, 2, 4, 2], [0, 7, 3, 3], [3, 3, 0, 2],
                      [1, 3, 2, 2]], 34),)

    known_lcs = (('AACCTTGG', 'ACACTGTGA', 'AACTGG'),)

    known_dag = ((0, 4, {0: {1: 7, 2: 4}, 1: {4: 1}, 2: {3: 2}, 3: {4: 3}},
                  9, [0, 2, 3, 4]),)

    known_align = (('PA', 'APA', 6), ('PLEASANTLY', 'MEANLY', 8))

    known_local = (('MEANLY', 'PENALTY', 15),)

    known_edit = (('MAP', 'CARP', 2), ('PLEASANTLY', 'MEANLY', 5))

    known_fitting = (('GTAGGCTTAAGGTTA', 'TAGATA', 2),)

    known_overlap = (('PAWHEAE', 'HEAGAWGHEE', 1),)

    known_affine = (('A', 'NA', -7), ('NA', 'A', -7),
                    ('PLAAN', 'PAN', 5), ('PAN', 'PLAAN', 5),
                    ('PRTEINS', 'PRTWPSEIN', 8))

    known_middle_edge = (('A', '', ((0, 0), 'v')),
                         ('PA', 'A', ((1, 0), 'd')),
                         ('PA', 'AA', ((1, 1), 'd')),
                         ('A', 'AP', ((1, 1), 'h')),
                         ('PAAA', 'APA', ((2, 1), 'd')),
                         ('PLEASANTLY', 'MEASNLY', ((4, 3), 'd')))

    def test_coins(self):
        """Test min-change finder"""
        for value, denoms, coins in self.known_coins:
            result = dynamic_practice.min_coins(value, denoms)
            self.assertEqual(result, coins)

    def test_length(self):
        """Test max-length finder"""
        for down, right, length in self.known_length:
            result = dynamic_practice.max_length(down, right)
            self.assertEqual(result, length)

    def test_lcs(self):
        """Test LCS finder"""
        for one, two, lcs in self.known_lcs:
            result = dynamic_practice.lcs(one, two)
            self.assertEqual(len(result), len(lcs))

    def test_dag(self):
        """Test generalized DAG"""
        for source, sink, paths, length, path in self.known_dag:
            graph = dag.DAG(source, sink, paths)
            result_l, result_p = graph.longest_path()
            self.assertEqual(result_l, length)
            self.assertEqual(result_p, path)

    def test_align(self):
        """Test global aligner"""
        score_matrix = aligner.read_score_matrix('blossom.txt')
        for one, two, score in self.known_align:
            result_s, result_a = aligner.global_align(one, two, score_matrix)
            self.assertEqual(result_s, score)

    def test_local(self):
        """Test local aligner"""
        score_matrix = aligner.read_score_matrix('pam.txt')
        for one, two, score in self.known_local:
            result_s, result_a = aligner.local_align(one, two, score_matrix)
            self.assertEqual(result_s, score)

    def test_edit(self):
        """Test edit distance calculator"""
        for one, two, distance in self.known_edit:
            result = aligner.find_edit_distance(one, two)
            self.assertEqual(result, distance)

    def test_fitting(self):
        """Test fitting aligner"""
        for long, short, score in self.known_fitting:
            result_s, result_a = aligner.fitting_align(long, short)
            self.assertEqual(result_s, score)

    def test_overlap(self):
        """Test overlap aligner"""
        for before, after, score in self.known_overlap:
            result_s, result_a = aligner.overlap_align(before, after)
            self.assertEqual(result_s, score)

    def test_affine(self):
        """Test affine aligner"""
        score_matrix = aligner.read_score_matrix('blossom.txt')
        for one, two, score in self.known_affine:
            result_s, result_a = aligner.affine_align(one, two, -11, -1,
                                                      score_matrix)
            self.assertEqual(result_s, score)

    def test_middle_edge(self):
        """Test middle edge finder"""
        score_matrix = aligner.read_score_matrix('blossom.txt')
        for one, two, edge in self.known_middle_edge:
            print("testing ", one, " and ", two,
                  " which should have edge ", edge)
            result = aligner.find_middle_edge(one, two, score_matrix, -5)            
            self.assertEqual(result, edge)

if __name__ == '__main__':
    unittest.main()

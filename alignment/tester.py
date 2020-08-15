import dynamic_practice
import dag
import protein_aligner
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
        for one, two, score in self.known_align:
            score_matrix = protein_aligner.read_score_matrix('blossom.txt')
            aligner = protein_aligner.ProteinAligner(one, two, score_matrix, -5)
            result_s, result_a = aligner.align()
            self.assertEqual(result_s, score)

if __name__ == '__main__':
    unittest.main()

import assembler
import unittest

class Tester(unittest.TestCase):
    known_comp = (('CAATCCAAC', 5,
                   ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']),)

    known_comp_failure = (('', 5),
                          ('CAAG', 5),
                          ('ATTC', 0))

    known_path_to_DNA = ((['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'],
                         'ACCGAAGCT'),)

    known_pattern_failure = (['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', ''],
                                 [], ['', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'],
                                 ['ACCA', 'CCGAA', 'CGAAG', 'GAAGC'],
                                 ['ACCGA', 'CCGAA', 'CGAAG', 'GA'])

    known_overlap = ((['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT'],
                      {'AGGCA': ['GGCAT'], 'CATGC': ['ATGCG'],
                       'GCATG': ['CATGC'], 'GGCAT': ['GCATG']}),)

    known_self_overlap = ((['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG',
                            'GGAG'],
                           {'AGG': ['GGG'], 'CAG': ['AGG', 'AGG'],
                            'GAG': ['AGG'], 'GGA': ['GAG'],
                            'GGG': ['GGG', 'GGA']}),)

    def test_comp(self):
        """Composition of a string should be found & sorted correctly"""

        for DNA, pat_len, comp in self.known_comp:
            result = assembler.comp(DNA, pat_len)
            self.assertEqual(comp, result)

    def test_comp_failure(self):
        """Composition finder should error on bad input"""

        for DNA, pat_len in self.known_comp_failure:
            self.assertRaises(ValueError, assembler.comp, DNA, pat_len)

    def test_path_to_DNA(self):
        """DNA should be properly reconstructed from path"""

        for path, DNA in self.known_path_to_DNA:
            result = assembler.path_to_DNA(path)
            self.assertEqual(DNA, result)

    def test_path_to_DNA_failure(self):
        """DNA-from-path should error on bad input"""

        for path in self.known_pattern_failure:
            self.assertRaises(ValueError, assembler.path_to_DNA, path)

    def test_overlap(self):
        """Overlap graph should be properly constructed"""

        for pats, graph in self.known_overlap:
            result = assembler.overlap_graph(pats)
            self.assertEqual(graph, result)

    def test_overlap_failure(self):
        """Overlap graph constructor should error on bad input"""

        for pats in self.known_pattern_failure:
            self.assertRaises(ValueError, assembler.overlap_graph, pats)

    def test_self_overlap(self):
        """Self-overlap graph should be properly constructed"""

        for pats, graph in self.known_self_overlap:
            result = assembler.self_overlap_graph(pats)
            self.assertEqual(graph, result)

    def test_self_overlap_failure(self):
        """Self-overlap graph constructor should error on bad input"""

        for pats in self.known_pattern_failure:
            self.assertRaises(ValueError, assembler.self_overlap_graph, pats)

if __name__ == '__main__':
    unittest.main()

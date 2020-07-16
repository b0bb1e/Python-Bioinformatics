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

    known_graph_to_cycle = (({'0': ['3'], '1': ['0'], '2': ['1', '6'],
                              '3': ['2'], '4': ['2'], '5': ['4'],
                              '6': ['5', '8'], '7': ['9'], '8': ['7'],
                              '9': ['6']},
                             ['3', '2', '6', '8', '7', '9', '6', '5', '4', '2',
                              '1', '0', '3']),)

    known_graph_to_path = (({'0': ['2'], '1': ['3'], '2': ['1'],
                             '3': ['0', '4'], '6': ['3', '7'], '7': ['8'],
                             '8': ['9'], '9': ['6']},
                            ['6', '7', '8', '9', '6', '3', '0', '2', '1', '3',
                             '4']),)

    known_assemble = ((['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC'],
                       'GGCTTACCA'),)

    known_assemble_read_pairs = (([('GAGA', 'TTGA'), ('TCGT', 'GATG'),
                                   ('CGTG', 'ATGT'), ('TGGT', 'TGAG'),
                                   ('GTGA', 'TGTT'), ('GTGG', 'GTGA'),
                                   ('TGAG', 'GTTG'), ('GGTC', 'GAGA'),
                                   ('GTCG', 'AGAT')], 2, 'GTGGTCGTGAGATGTTGA'),)
    
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

    def test_graph_to_cycle(self):
        """Eulerian cycle maker should construct cycles properly"""

        for graph, cycle in self.known_graph_to_cycle:
            result = assembler.graph_to_cycle(graph)
            self.assertEqual(cycle, result)

    def test_graph_to_path(self):
        """Eulerian path maker should construct paths properly"""


        for graph, path in self.known_graph_to_path:
            result = assembler.graph_to_path(graph)
            self.assertEqual(path, result)

    def test_assemble(self):
        """DNA assembler should assemble correctly"""

        for pats, DNA in self.known_assemble:
            result = assembler.assemble(pats)
            self.assertEqual(DNA, result)

    def test_assemble_failure(self):
        """DNA assembler should error on bad input"""

        for pats in self.known_pattern_failure:
            self.assertRaises(ValueError, assembler.assemble, pats)

    def test_assemble_read_pairs(self):
        """DNA read-pair assembler should assemble correctly"""

        for pairs, dist, DNA in self.known_assemble_read_pairs:
            result = assembler.assemble_read_pairs(pairs, dist)
            self.assertEqual(DNA, result)

if __name__ == '__main__':
    unittest.main()

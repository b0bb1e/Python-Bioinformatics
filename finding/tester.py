import freq_finder
import unittest

class Tester(unittest.TestCase):
    """Unit tester for freq_finder.py"""

    # test cases for find_starts
    known_starts = (('AGGGTCAGCGATCA', 'TCA', [4, 11]),
                    ('AAAAAAAAA', 'AA', [0, 1, 2, 3, 4, 5, 6, 7]),
                    ('GAGCTGCACCCC', 'AAA', []),
                    ('CATAT', 'CATATAT', []))

    # test cases for calc_freq
    known_freq = (('ACGAGTAC', 2, [0, 2, 1, 0, 0, 0, 1, 0,
                                   1, 0, 0, 1, 1, 0, 0, 0]),
                  ('AAAAAA', 1, [6, 0, 0, 0]))

    # test cases for find_most_freq
    known_most_freq = (('AGACTCAGCTTAG', 2, 2, True, ['AG', 'CT']),
                       ('AGACTCAGCTTAG', 2, 2, False, ['AG']),
                       ('AGACTCAGCTTAG', 2, 3, False, ['AG']),
                       ('AGACTCAGCTTAG', 2, 4, True, []))

    # test cases for find_clumps
    known_clumps = (('AGCATGATGTGACAGTAC', 2, 3, 8, ['TG']),
                    ('ATAATGCTGTGACATTAT', 2, 3, 8, ['TG']))

    # test cases for rev_comp
    known_rev_comp = (('AAAACCCGGT', 'ACCGGGTTTT'),
                      ('ACACAC', 'GTGTGT'))

    # test cases for find_approx_starts
    known_approx_starts = (('CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCAC'
                            + 'GGTACGGACGTCAATCAAAT', 'ATTCTGGA', 3,
                            [6, 7, 26, 27]),
                           ('TTTTTTAAATTTTAAATTTTTT', 'AAA', 2,
                            [4, 5, 6, 7, 8, 11, 12, 13, 14, 15]),
                           ('GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATT'
                            + 'GGCGATGAAACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCC'
                            + 'GCTATCGCATTTTGATGCGCGCCGCGTCGATT', 'GAGCGCTGG',
                            2, [0, 30, 66]),
                           ('CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTC'
                            + 'TGCATATAAGTAGTGAAGGTATAGAAACCCGTTCAAGCCCGCAGCGG'
                            + 'TAAAACCGAGAACCATGATGAATGCACGGCGATTGCGCCATAATCCA'
                            + 'AACA', 'AATCCTTTCA', 3, [3, 36, 74, 137]),
                           ('CCGTCATCCGTCATCCTCGCCACGTTGGCATGCATTCCGTCATCCCGTC'
                            + 'AGGCATACTTCTGCATATAAGTACAAACATCCGTCATGTCAAAGGGA'
                            + 'GCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGC'
                            , 'CCGTCATCC', 3, [0, 7, 36, 44, 48, 72, 79, 112]),
                           ('AAAAAA', 'TTT', 3, [0, 1, 2, 3]),
                           ('CCACCT', 'CCA', 0, [0]))

    # test cases for find_most_approx_freq
    known_approx_freq = (('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1,
                          ['ATGC','ATGT', 'GATG']),
                         ('AAAAAAAAAA', 2, 1,
                          ['AA', 'AC', 'AG', 'CA', 'AT', 'GA', 'TA']),
                         ('AGTCAGTC', 4, 2,
                          ['TCTC', 'CGGC', 'AAGC', 'TGTG', 'GGCC', 'AGGT',
                           'ATCC', 'ACTG', 'ACAC', 'AGAG', 'ATTA', 'TGAC',
                           'AATT', 'CGTT', 'GTTC', 'GGTA', 'AGCA', 'CATC']),
                         ('AATTAATTGGTAGGTAGGTA', 4, 0, ['GGTA']))

    def test_num_to_DNA_to_num(self):
        """number -> DNA string conversion should be reversible"""

        for chars in range(1, 8):
            # loop over all VALID numbers for this length
            for i in range(4 ** chars):
                DNA = freq_finder.num_to_DNA(i, chars)
                num = freq_finder.DNA_to_num(DNA)
                self.assertEqual(i, num)
                new_DNA = freq_finder.num_to_DNA(num, chars)
                self.assertEqual(DNA, new_DNA)

    def test_num_to_DNA_failure(self):
        """number -> DNA string conversion should error on bad input"""

        self.assertRaises(ValueError, freq_finder.num_to_DNA, -1, 10)
        self.assertRaises(ValueError, freq_finder.num_to_DNA, 17, 2)
        self.assertRaises(ValueError, freq_finder.num_to_DNA, 3, 0)

    def test_DNA_to_num_failure(self):
        """DNA string -> number conversion should error on bad input"""

        self.assertRaises(ValueError, freq_finder.DNA_to_num, ' ')
        self.assertRaises(ValueError, freq_finder.DNA_to_num, '')
        self.assertRaises(ValueError, freq_finder.DNA_to_num, 'GAGTAB')
        self.assertRaises(ValueError, freq_finder.DNA_to_num, 'BGAGCA')

    def test_start_finder(self):
        """Start indicies of a substring should be collected correctly"""

        for DNA, pattern, starts in self.known_starts:
            result = freq_finder.find_starts(DNA, pattern)
            self.assertEqual(starts, result)

    def test_start_finder_failure(self):
        """Start-finder should error on bad input"""

        self.assertRaises(ValueError, freq_finder.find_starts, '', '')
        self.assertRaises(ValueError, freq_finder.find_starts, 'A', '')
        self.assertRaises(ValueError, freq_finder.find_starts, '', 'A')

    def test_freq_finder(self):
        """Frequency-finder should calculate frequencies correctly"""

        for DNA, pat_len, freq in self.known_freq:
            result = freq_finder.calc_freq(DNA, pat_len)
            self.assertEqual(freq, result)

    def test_freq_finder_failure(self):
        """Frequency-finder should error on bad input"""

        self.assertRaises(ValueError, freq_finder.calc_freq, '', 4)
        self.assertRaises(ValueError, freq_finder.calc_freq, '', -4)
        self.assertRaises(ValueError, freq_finder.calc_freq, 'AGTA', 0)
        self.assertRaises(ValueError, freq_finder.calc_freq, 'AGYTA', 2)
        self.assertRaises(ValueError, freq_finder.calc_freq, ' AGTA', 0)

    def test_freq_pattern_finder(self):
        """Frequent-pattern-finder should find patterns correctly"""

        for DNA, pat_len, min_times, all_above, freq in self.known_most_freq:
            result = freq_finder.find_most_freq(DNA, pat_len, min_times,
                                                all_above)
            self.assertEqual(freq, result)

    def test_freq_pattern_finder_failure(self):
        """Frequent-pattern-finder should error on bad input"""
        
        self.assertRaises(ValueError, freq_finder.find_most_freq, '', 3)
        self.assertRaises(ValueError, freq_finder.find_most_freq, 'AGTAG', 0)
        self.assertRaises(ValueError, freq_finder.find_most_freq, 'AGAT', 3, 1)
        self.assertRaises(ValueError, freq_finder.find_most_freq, 'JAGTAG', 2)
        self.assertRaises(ValueError, freq_finder.find_most_freq, 'AGTWAG', 0)

    def test_clump_finder(self):
        """Clump-finder should find clumped patterns correctly"""

        for DNA, pat_len, min_times, window_len, clumps in self.known_clumps:
            result = freq_finder.find_clumps(DNA, pat_len, min_times,
                                             window_len)
            self.assertEqual(clumps, result)

    def test_clump_finder_failure(self):
        """Clump-finder should error on bad input"""

        self.assertRaises(ValueError, freq_finder.find_clumps, '', 3, 2, 10)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'AGGA', 0, 2, 3)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'ACA', 3, 1, 10)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'ACAT', 3, 2, 6)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'ACAT', 2, 2, 2)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'OACAT', 2, 2, 4)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'ACHAT', 2, 2, 4)
        
    def test_rev_comp(self):
        """Reverse complement should be found correctly"""
        
        for pat, rev in self.known_rev_comp:
            result = freq_finder.rev_comp(pat)
            self.assertEqual(rev, result)

    def test_rev_comp_failure(self):
        """Reverse complement maker should error on bad input"""

        self.assertRaises(ValueError, freq_finder.rev_comp, '')
        self.assertRaises(ValueError, freq_finder.rev_comp, 'L')
        self.assertRaises(ValueError, freq_finder.rev_comp, 'ACGTR')

    def test_approx_start_finder(self):
        """Approximate starting positions should be found correctly"""

        for DNA, pat, dist, starts in self.known_approx_starts:
            result = freq_finder.find_approx_starts(DNA, pat, dist)
            self.assertEqual(starts, result)

    def test_approx_start_finder_failure(self):
        """Approximate starting position finder should error on bad input"""

        self.assertRaises(ValueError, freq_finder.find_approx_starts,
                         'AAAAAA', 'TTT', -1)
        self.assertRaises(ValueError, freq_finder.find_approx_starts,
                         '', 'TTT', 1)
        self.assertRaises(ValueError, freq_finder.find_approx_starts,
                         'AAAAA', '', 1)

    def test_approx_freq_pattern_finder(self):
        """Most-frequent approximate patterns should be found correctly"""

        for DNA, pat_len, dist, freq in self.known_approx_freq:
            result = freq_finder.find_most_approx_freq(DNA, pat_len, dist)
            self.assertCountEqual(freq, result)

    def test_approx_freq_pattern_finder_failure(self):
        """Most-frequent approximate pattern finder should error on bad input"""

        self.assertRaises(ValueError, freq_finder.find_most_approx_freq,
                          'AATTAATTGGTAGGTAGGTA', 4, -1)
        self.assertRaises(ValueError, freq_finder.find_most_approx_freq,
                          'AATTAATTGGTAGGTAGGTA', 0, 0)
        self.assertRaises(ValueError, freq_finder.find_most_approx_freq,
                          'AATTAATYTGGTAGGTAGGTA', 4, 2)
        self.assertRaises(ValueError, freq_finder.find_most_approx_freq,
                          'EAATTAATTGGTAGGTAGGTA', 4, 1)
        self.assertRaises(ValueError, freq_finder.find_most_approx_freq,
                          'ATA', 4, 1)
                          
if __name__ == '__main__':
    unittest.main()

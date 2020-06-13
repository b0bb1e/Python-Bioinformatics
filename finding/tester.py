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

    # test cases for find_freq
    known_most_freq = (('AGACTCAGCTTAG', 2, 2, True, ['AG', 'CT']),
                       ('AGACTCAGCTTAG', 2, 2, False, ['AG']),
                       ('AGACTCAGCTTAG', 2, 3, False, ['AG']),
                       ('AGACTCAGCTTAG', 2, 4, True, []))

    # test cases for find_clumps
    known_clumps = (('AGCATGATGTGACAGTAC', 2, 3, 8, ['TG']),
                    ('ATAATGCTGTGACATTAT', 2, 3, 8, ['TG']))

    def test_num_to_DNA_to_num(self):
        """number -> DNA string conversion should be reversible"""

        for chars in range(8):
            # loop over all VALID numbers for this length
            for i in range(4 ** chars):
                DNA = freq_finder.num_to_DNA(i, chars)
                num = freq_finder.DNA_to_num(DNA)
                self.assertEqual(i, num)
                new_DNA = freq_finder.num_to_DNA(num, chars)
                self.assertEqual(DNA, new_DNA)

    def test_start_finder(self):
        """Start indicies of a substring should be collected correctly"""

        for DNA, pattern, starts in self.known_starts:
            result = freq_finder.find_starts(DNA, pattern)
            self.assertEqual(starts, result)

    def test_start_finder_failure(self):
        """Start-finder should error on empty strings"""

        self.assertRaises(ValueError, freq_finder.find_starts, '', '')
        self.assertRaises(ValueError, freq_finder.find_starts, 'A', '')
        self.assertRaises(ValueError, freq_finder.find_starts, '', 'A')

    def test_freq_finder(self):
        """Frequency-finder should calculate frequencies correctly"""

        for DNA, pat_len, freq in self.known_freq:
            result = freq_finder.calc_freq(DNA, pat_len)
            self.assertEqual(freq, result)

    def test_freq_finder_failure(self):
        """Frequency-finder should error on bad inputs"""

        self.assertRaises(ValueError, freq_finder.calc_freq, '', 4)
        self.assertRaises(ValueError, freq_finder.calc_freq, '', -4)
        self.assertRaises(ValueError, freq_finder.calc_freq, 'AGTA', 0)

    def test_freq_pattern_finder(self):
        """Frequent-pattern-finder should find
        the most-frequent patterns correctly
        """

        for DNA, pat_len, min_times, all_above, freq in self.known_most_freq:
            result = freq_finder.find_freq(DNA, pat_len, min_times, all_above)
            self.assertEqual(freq, result)

    def test_freq_pattern_finder_failure(self):
        """Frequent-pattern-finder should error on bad input"""
        
        self.assertRaises(ValueError, freq_finder.find_freq, '', 3)
        self.assertRaises(ValueError, freq_finder.find_freq, 'AGTAG', 0)
        self.assertRaises(ValueError, freq_finder.find_freq, 'AGATA', 3, 1)

    def test_clump_finder(self):
        """Clump-finder should find clumped patterns correctly"""

        for DNA, pat_len, min_times, window_length, clumps in self.known_clumps:
            result = freq_finder.find_clumps(DNA, pat_len, min_times,
                                             window_length)
            self.assertEqual(clumps, result)

    def test_clump_finder_failure(self):
        """Clump-finder should error on bad input"""

        self.assertRaises(ValueError, freq_finder.find_clumps, '', 3, 2, 10)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'AGGAG', 0, 2, 10)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'ACA', 3, 1, 10)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'ACAT', 3, 2, 6)
        self.assertRaises(ValueError, freq_finder.find_clumps, 'ACAT', 2, 2, 2)
        
if __name__ == '__main__':
    unittest.main()

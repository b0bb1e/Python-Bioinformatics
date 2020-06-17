import motifs
import unittest

class Tester(unittest.TestCase):
    """Unit tester for motifs.py"""

    # test cases for brute_finder
    known_brute = ((['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'], 3, 1,
                    {'ATA', 'ATT', 'GTT', 'TTT'}),
                   (['ACGT', 'ACGT', 'ACGT'], 3, 0, {'ACG', 'CGT'}),
                   (['AAAAA', 'AAAAA', 'AAAAA'], 3, 1,
                    {'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'AGA', 'ATA', 'CAA',
                     'GAA', 'TAA'}),
                   (['AAAAA', 'AAAAA', 'AACAA'], 3, 0, set()),
                   (['AACAA', 'AAAAA', 'AACAA'], 3, 0, set()))

    # test cases for median_string
    known_median = ((['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG',
                      'GCTGAGCACCGG', 'AGTACGGGACAG'], 3, {'ACG', 'GAC'}),
                    (['ACGT', 'ACGT', 'ACGT'], 3, {'ACG', 'CGT'}),
                    (['ATA', 'ACA', 'AGA', 'AAT', 'AAC'], 3, {'AAA'}))

    def test_brute(self):
        """Brute-force motif finder should find exact motifs"""

        for DNAs, pat_len, dist, correct in self.known_brute:
            result = motifs.brute_finder(DNAs, pat_len, dist)
            self.assertEqual(correct, result)

    def test_brute_failure(self):
        """Brute-force motif finder should error on bad input"""
        
        self.assertRaises(ValueError, motifs.brute_finder,
                          ['ACGT', 'ACGR', 'ACGT'], 3, 0)
        self.assertRaises(ValueError, motifs.brute_finder,
                          ['YACG', 'ACGT', 'ACGT'], 3, 0)
        self.assertRaises(ValueError, motifs.brute_finder,
                          ['ACGT', 'ACGT', 'ACGT'], 3, -1)
        self.assertRaises(ValueError, motifs.brute_finder,
                          ['ACGT', 'ACGT', 'ACGT'], 0, 0)
        self.assertRaises(ValueError, motifs.brute_finder,
                          ['ACGT', '', 'ACGT'], 3, 0)
        self.assertRaises(ValueError, motifs.brute_finder,
                          [''], 3, 0)
        self.assertRaises(ValueError, motifs.brute_finder,
                          ['ACGT'], 3, 0)
        self.assertRaises(ValueError, motifs.brute_finder,
                          [], 3, 0)

    def test_median(self):
        """Median string finder should find a valid best string"""

        for DNAs, pat_len, corrects in self.known_median:
            result = motifs.median_string(DNAs, pat_len)
            self.assertIn(result, corrects)

    def test_median_failure(self):
        """Median string finder should error on bad input"""

        self.assertRaises(ValueError, motifs.median_string,
                          ['ACGT', 'ACGT', 'ACGT'], 0)
        self.assertRaises(ValueError, motifs.median_string, [], 3)
        self.assertRaises(ValueError, motifs.median_string,
                          ['ACGT', 'ACGT', ''], 3)
        self.assertRaises(ValueError, motifs.median_string,
                          ['ACGT'], 3)

if __name__ == '__main__':
    unittest.main()

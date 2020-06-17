import motifs
import unittest

class Tester(unittest.TestCase):
    known_brute = ((['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'], 3, 1,
                    {'ATA', 'ATT', 'GTT', 'TTT'}),
                   (['ACGT', 'ACGT', 'ACGT'], 3, 0, {'ACG', 'CGT'}),
                   (['AAAAA', 'AAAAA', 'AAAAA'], 3, 1,
                    {'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'AGA', 'ATA', 'CAA',
                     'GAA', 'TAA'}),
                   (['AAAAA', 'AAAAA', 'AACAA'], 3, 0, set()),
                   (['AACAA', 'AAAAA', 'AAaAA'], 3, 0, set()))

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

if __name__ == '__main__':
    unittest.main()

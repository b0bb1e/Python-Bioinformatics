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

    known_brute_failure = ((['ACGT', 'ACGR', 'ACGT'], 3, 0),
                           (['YACG', 'ACGT', 'ACGT'], 3, 0),
                           (['ACGT', 'ACGT', 'ACGT'], 3, -1),
                           (['ACGT', 'ACGT', 'ACGT'], 0, 0),
                           (['ACGT', '', 'ACGT'], 3, 0),
                           ([''], 3, 0),
                           (['ACGT'], 3, 0),
                           ([], 3, 0))

    # test cases for median_string
    known_median = ((['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG',
                      'GCTGAGCACCGG', 'AGTACGGGACAG'], 3, {'ACG', 'GAC'}),
                    (['ACGT', 'ACGT', 'ACGT'], 3, {'ACG', 'CGT'}),
                    (['ATA', 'ACA', 'AGA', 'AAT', 'AAC'], 3, {'AAA'}))

    known_median_failure = ((['ACGT', 'ACGT', 'ACGT'], 0),
                            ([], 3),
                            (['ACGT', 'ACGT', ''], 3),
                            (['ACGT'], 3))

    # test cases for greedy_finder
    known_greedy = ((['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC',
                      'CACGTCAATCAC', 'CAATAATATTCG'], 3,
                     ['TTC', 'ATC', 'TTC', 'ATC', 'TTC']),
                    (['AGGCGGCACATCATTATCGATAACGATTCGCCGCATTGCC',
                      'ATCCGTCATCGAATAACTGACACCTGCTCTGGCACCGCTC',
                      'AAGCGTCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT',
                      'AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
                      'AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT',
                      'AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT',
                      'AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG',
                      'AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA'], 5,
                     ['AGGCG', 'ATCCG', 'AAGCG', 'AGTCG', 'AACCG', 'AGGCG',
                      'AGGCG', 'AGGCG']),
                    (['GCACATCATTAAACGATTCGCCGCATTGCCTCGATAGGCG',
                      'TCATAACTGACACCTGCTCTGGCACCGCTCATCCGTCGAA',
                      'AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTCGGC',
                      'AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
                      'AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT',
                      'AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT',
                      'AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTGAGGCG',
                      'AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCG'], 5,
                     ['AGGCG', 'TGGCA', 'AAGCG', 'AGGCA', 'CGGCA', 'AGGCG',
                      'AGGCG', 'AGGCG']),
                    (['GCACATCATTATCGATAACGATTCATTGCCAGGCGGCCGC',
                      'TCATCGAATAACTGACACCTGCTCTGGCTCATCCGACCGC',
                      'TCGGCGGTATAGCCAGATAGTGCCAATAATTTCCTAAGCG',
                      'GTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTGAGTCG',
                      'GACGGCAACTACGGTTACAACGCAGCAAGAATATTAACCG',
                      'TCTGTTGTTGCTAACACCGTTAAGCGACGGCAACTAGGCG',
                      'GCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTGAAGCG',
                      'AAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAAAATTG'], 5,
                     ['GGCGG', 'GGCTC', 'GGCGG', 'GGCAG', 'GACGG', 'GACGG',
                      'GGCGC', 'GGCGC']))

    # failure cases for greedy or random
    known_g_or_r_failure = ((['GCCYAA', 'GGCCTG', 'AACCTA', 'TTCCTT'], 3),
                            (['GCCCAA', 'RGCCTG', 'AACCTA', 'TTCCTT'], 3),
                            (['GCCCAA', 'GGCCTG', 'AACCTA', 'TTCCTT'], 0),
                            (['GCCCAA', 'GGCCTG', '', 'TTCCTT'], 3),
                            (['GCCCAA'], 3),
                            ([''], 3),
                            ([], 3))

    # test cases for random_finder
    known_random = ((['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                      'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                      'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                      'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                      'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'], 8,
                     ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG',
                      'TCCACGTG']),
                    (['AATTGGCACATCATTATCGATAACGATTCGCCGCATTGCC',
                      'GGTTAACATCGAATAACTGACACCTGCTCTGGCACCGCTC',
                      'AATTGGCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT',
                      'GGTTAATGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
                      'AATTGGACGGCAACTACGGTTACAACGCAGCAAGAATATT',
                      'GGTTAACTGTTGTTGCTAACACCGTTAAGCGACGGCAACT',
                      'AATTGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG',
                      'GGTTAAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA'], 6,
                     ['CGATAA', 'GGTTAA', 'GGTATA', 'GGTTAA', 'GGTTAC',
                      'GGTTAA', 'GGCCAA', 'GGTTAA']),
                    (['GCACATCATTAAACGATTCGCCGCATTGCCTCGATTAACC',
                      'TCATAACTGACACCTGCTCTGGCACCGCTCATCCAAGGCC',
                      'AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTAACC',
                      'AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAAGGCC',
                      'AACCGGACGGCAACTACGGTTACAACGCAGCAAGTTAACC',
                      'AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGAAGGCC',
                      'AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTTTAACC',
                      'AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCC'], 6,
                     ['TTAACC', 'ATAACT', 'TTAACC', 'TGAAGT', 'TTAACC',
                      'TTAAGC', 'TTAACC', 'TGAACA']))

    def test_brute(self):
        """Brute-force motif finder should find exact motifs"""

        test_num = 1
        for DNAs, pat_len, dist, correct in self.known_brute:
            result = motifs.brute_finder(DNAs, pat_len, dist)
            self.assertEqual(correct, result, msg=('Test #' + str(test_num)))
            test_num += 1

    def test_brute_failure(self):
        """Brute-force motif finder should error on bad input"""

        for DNAs, pat_len, dist in self.known_brute_failure:
            self.assertRaises(ValueError, motifs.brute_finder,
                              DNAs, pat_len, dist)

    def test_median(self):
        """Median string finder should find a valid best string"""

        for DNAs, pat_len, corrects in self.known_median:
            result = motifs.median_string(DNAs, pat_len)
            self.assertIn(result, corrects)

    def test_median_failure(self):
        """Median string finder should error on bad input"""

        for DNAs, pat_len in self.known_median_failure:
            self.assertRaises(ValueError, motifs.median_string, DNAs, pat_len)

    def test_greedy(self):
        """Greedy motif finder should find expected motifs"""

        test_num = 1
        for DNAs, pat_len, correct in self.known_greedy:
            result = motifs.greedy_finder(DNAs, pat_len)
            self.assertEqual(correct, result, msg=('Test #' + str(test_num)))
            test_num += 1

    def test_greedy_failure(self):
        """Greedy motif finder should error on bad input"""
        
        for DNAs, pat_len in self.known_g_or_r_failure:
            self.assertRaises(ValueError, motifs.greedy_finder, DNAs, pat_len)
#Set up before cutting (7) -> PREPARE
    def test_random(self):
        """Random motif finder should find expected motifs"""

        test_num = 1
        for DNAs, pat_len, correct in self.known_random:
            result = motifs.random_finder(DNAs, pat_len)
            self.assertEqual(correct, result, msg=('Test #' + str(test_num)))
            test_num += 1

    def test_random_failure(self):
        """Random motif finder should error on bad input"""
        
        for DNAs, pat_len in self.known_g_or_r_failure:
            self.assertRaises(ValueError, motifs.random_finder, DNAs, pat_len)

    def test_sampler(self):
        """Random-sampler motif finder should find expected motifs"""

        test_num = 1
        for DNAs, pat_len, correct in self.known_random:
            result = motifs.sampler_finder(DNAs, pat_len)
            self.assertEqual(correct, result, msg=('Test #' + str(test_num)))
            test_num += 1

    def test_sampler_failure(self):
        """Random-sampler motif finder should error on bad input"""
        
        for DNAs, pat_len in self.known_g_or_r_failure:
            self.assertRaises(ValueError, motifs.sampler_finder, DNAs, pat_len)
            
if __name__ == '__main__':
    unittest.main()

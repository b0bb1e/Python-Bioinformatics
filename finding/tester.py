import frequency_finder
import unittest

class Tester(unittest.TestCase):
    num_to_DNA = ((1054373062, 'TGTCGACTCATACG'),)
      #            (2044610837, 'TGCTCTGCAGCACCC'),
      #            (1600148523, 'CTTCGAACCCAAGGT'),
      #            (26552272, 'GCCCAGCTTCAA'),
      #            (1976879007, 'TCCTCCATACTGCTT'),
      #            (777314701, 'TGCCCATGACGATC'),
      #            (58956729, 'GAATGCGTGTGC'),
      #            (602462345, 'ATTGGATCGAGAGC'),
      #            (100903030, 'GAAATGGGACTCG'),
      #            (1276043580, 'ATAAATGTGACATTA'),
      #            (775604526, 'TGATGGTAGCAGTG'),
      #            (811464132, 'AACCTCTTCCTACA'),
      #            (2110241176, 'TTCTACTGTGCGCGA'),
      #            (1175370808, 'ACGAATGGTTAATGA'),
      #            (1479931973, 'CGAATCCTTGACACC'),
      #            (1742589314, 'GCTTCTCTATCGAAG'),
      #            (711061844, 'GGCGACTTACCCCA'),
      #            (1485624572, 'CGAGATATCCATTTA'),
      #           (995480382, 'GTCCCCTCATATTG'),
      #            (1250341674, 'AGGGACGGTATAGGG'),
      #            (583765934, 'AGTAGTGATTGGTG'),
      #            (1476132109, 'CCTTTGTTTTCAATC'),
      #            (8825740, 'ACGGGGTGATA'),
      #            (1301300237, 'ATCGCAACACAAATC'),
      #            (1694539068, 'GCCAAAAGCTCATTA'))

    def test_num_to_DNA(self):
        '''numbers should be converted into DNA strings correctly'''
        for num, DNA in self.num_to_DNA:
            result = frequency_finder.num_to_DNA(num, len(DNA))
            self.assertEqual(DNA, result)
            pat = frequency_finder.DNA_to_num(DNA)
            result = frequency_finder.num_to_DNA(pat, len(DNA))
            

    #def test_DNA_to_num(self):
    #    '''DNA strings should be converted into numbers correctly'''
    #    for num, DNA in self.num_to_DNA:
    #        result = frequency_finder.DNA_to_num(DNA)
    #        self.assertEqual(num, result)

    #def test_num_to_DNA_to_num(self):
    #    '''number -> DNA string conversion should be reversible'''
    #    for i in range(4 ** 7):
            #print('i=' + str(i))
    #        DNA = frequency_finder.num_to_DNA(i, 7)
            #print(DNA)
    #        num = frequency_finder.DNA_to_num(DNA)
            #print(num)
    #        self.assertEqual(i, num)
            
if __name__ == '__main__':
    unittest.main()

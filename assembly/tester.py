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

    known_path_to_DNA_failure = (['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', ''],
                                 [],
                                 ['', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT'],
                                 ['ACCA', 'CCGAA', 'CGAAG', 'GAAGC'],
                                 ['ACCGA', 'CCGAA', 'CGAAG', 'GA'])

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

        for path in self.known_path_to_DNA_failure:
            self.assertRaises(ValueError, assembler.path_to_DNA, path)

if __name__ == '__main__':
    unittest.main()

import sequencer
import spectrum
import unittest

class Tester(unittest.TestCase):
    """Tester for sequencer.py & spectrum.py"""
    
    known_RNA = (('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA',
                  'MAMAPRTEINSTRING'),)

    known_hidden_proteins = (('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGG'
                             + 'TGA', 'MA', ['ATGGCC', 'GGCCAT', 'ATGGCC']),)

    known_ideal_spectrum = (('LEQN', [0, 113, 114, 128, 129, 227, 242, 242, 257,
                                      355, 356, 370, 371, 484,]),)

    known_mass_count = ((1024, 14712706211),)

    known_sequencing = (([0, 113, 128, 186, 241, 299, 314, 427],
                         [[186, 128, 113], [186, 113, 128], [128, 186, 113],
                          [128, 113, 186], [113, 186, 128], [113, 128, 186]]),)

    known_score = (('NQEL', [0, 99, 113, 114, 128, 227, 257, 299, 355, 356,
                             370, 371, 484], 11),)

    known_leaderboard = (([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347,
                           389, 460,], 10, [113, 147, 71, 129]),)

    known_convolution = (([57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228,
                           265, 285, 299, 307, 323, 356, 364, 394, 422, 493],
                          60, 20, [99, 71, 137, 57, 72, 57]),)

    def test_RNA_sequencer(self):
        """RNA sequencer should work correctly"""
    
        for RNA, protein in self.known_RNA:
            result = sequencer.translate_RNA(RNA)
            self.assertEqual(result, protein)

    def test_protein_finder(self):
        """Encoding-finder should work correctly"""
    
        for DNA, protein, subs in self.known_hidden_proteins:
            result = sequencer.find_hidden_proteins(DNA, protein, True)
            self.assertCountEqual(result, subs)

    def test_ideal_spectrum(self):
        """Ideal spectrum calculator should work correctly"""

        for peptide, spect in self.known_ideal_spectrum:
            result = spectrum.ideal_spectrum_amino(peptide, True)
            self.assertEqual(result, spect)

    def test_mass_counting(self):
        """Peptides-with-mass counter should work correctly"""

        for mass, count in self.known_mass_count:
            result = spectrum.count_peptides_with_mass(mass)
            self.assertEqual(result, count)

    def test_sequencing(self):
        """Specturm sequencer should work correctly"""

        for spect, peptides in self.known_sequencing:
            result = spectrum.sequence(spect)
            self.assertEqual(result, peptides)

    def test_score(self):
        """Peptide scoreer should work correctly"""

        for peptide, spect, score in self.known_score:
            result = spectrum.score_amino(peptide, spect, True)
            self.assertEqual(result, score)

    def test_leaderboard(self):
        """Leaderboard sequencer should work correctly"""

        for spect, keep, peptide in self.known_leaderboard:
            result = spectrum.leaderboard_sequence(spect, keep)
            # checks that the lists are cyclic-equal
            self.assertIn(' '.join(map(str, result)),
                          ' '.join(map(str, peptide * 2)))

    def test_convolution(self):
        """Convolution sequencer should work correctly"""

        for spect, keep, num_acids, peptide in self.known_convolution:
            result = spectrum.convolution_sequence(spect, keep, num_acids)
            # checks that the lists are cyclic-equal
            self.assertIn(' '.join(map(str, result)),
                          ' '.join(map(str, peptide * 2)))

if __name__ == '__main__':
    unittest.main()
                

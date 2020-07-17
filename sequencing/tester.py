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
            self.assertEqual(result, spect)

if __name__ == '__main__':
    unittest.main()
                

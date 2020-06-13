import frequency_finder
import unittest

class Tester(unittest.TestCase):
    """Unit tester for frequency_finder.py"""
    
    def test_num_to_DNA_to_num(self):
        """number -> DNA string conversion should be reversible"""

        for chars in range(15):
            # loop over all VALID numbers for this length
            for i in range(4 ** chars):
                DNA = frequency_finder.num_to_DNA(i, chars)
                num = frequency_finder.DNA_to_num(DNA)
                self.assertEqual(i, num)
                new_DNA = frequency_finder.num_to_DNA(num, chars)
                self.assertEqual(DNA, new_DNA)
            
if __name__ == '__main__':
    unittest.main()

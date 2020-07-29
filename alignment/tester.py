import dynamic_practice
import unittest

class Tester(unittest.TestCase):
    known_coins = ((40, [1, 5, 10, 20, 25, 50], 2),)

    known_length = (([[1, 0, 2, 4, 3], [4, 6, 5, 2, 1], [4, 4, 5, 2, 1],
                      [5, 6, 8, 5, 3]],
                     [[3, 2, 4, 0], [3, 2, 4, 2], [0, 7, 3, 3], [3, 3, 0, 2],
                      [1, 3, 2, 2]], 34),)

    known_lcs = (('AACCTTGG', 'ACACTGTGA', 'AACTGG'),)

    def test_coins(self):
        for value, denoms, coins in self.known_coins:
            result = dynamic_practice.min_coins(value, denoms)
            self.assertEqual(result, coins)

    def test_length(self):
        for down, right, length in self.known_length:
            result = dynamic_practice.max_length(down, right)
            self.assertEqual(result, length)

    def test_lcs(self):
        for one, two, lcs in self.known_lcs:
            result = dynamic_practice.lcs(one, two)
            self.assertEqual(result, two)

if __name__ == '__main__':
    unittest.main()
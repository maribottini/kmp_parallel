import unittest
from kmp_parallel import create_kmp_table, kmp_search

class TestKMP(unittest.TestCase):

    def test_create_kmp_table(self):
        pattern = "BABABABAC"
        expected_table = [0, 0, 1, 2, 0, 1, 2, 3, 4]
        self.assertEqual(create_kmp_table(pattern).tolist(), expected_table)

    def test_kmp_search(self):
        genome = "ABABDABACDABABABABABACBCABAB"
        pattern = "BABABABAC"
        lock = None
        results = []
        kmp_search((genome, pattern, lock, results))
        expected_results = [(pattern, 10)]
        self.assertEqual(results, expected_results)

if __name__ == '__main__':
    unittest.main()

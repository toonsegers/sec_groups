import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import unittest

from mpyc.runtime import mpc
from demos.demo_basic_examples import suite1 


class CircuitSat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_test_suite1(self):
        verification = mpc.run(suite1())
        self.assertEqual(verification, True)


if __name__ == "__main__":
    unittest.main()

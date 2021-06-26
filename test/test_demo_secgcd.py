import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import unittest

from mpyc.runtime import mpc
from demos.demo_secgcd import (
    demo_secure_xgcd,
    demo_secure_division,
    demo_gcd_other,
    demo_secure_montgomery_exponentiation,
)


class CircuitSat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_demo_secure_xgcd(self):
        verification = mpc.run(demo_secure_xgcd(3))
        self.assertEqual(all(check == True for check in verification.values()), True)

    def test_demo_secure_division(self):
        verification = mpc.run(demo_secure_division(3))
        self.assertEqual(all(check == True for check in verification.values()), True)

    def test_demo_gcd_other(self):
        verification = mpc.run(demo_gcd_other(3))
        self.assertEqual(all(check == True for check in verification.values()), True)

    def test_demo_secure_montgomery_exponentiation(self):
        verification = mpc.run(demo_secure_montgomery_exponentiation(3))
        self.assertEqual(all(check == True for check in verification.values()), True)


if __name__ == "__main__":
    unittest.main()

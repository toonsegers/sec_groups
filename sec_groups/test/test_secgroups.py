import unittest
import os
import sys

project_root = sys.path.append(os.path.abspath("../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from mpyc.runtime import mpc
from tools.find_primes import find_safe_primes
import sec_groups.ellcurves as ell
import sec_groups.fingroups as fg
import sec_groups.secgroups as sfg
import sec_groups.pairing as pairing


class SecFiniteGrps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_sec_curves(self):
        curves = [
            fg.EllipticCurve(ell.ED25519, ell.ED_AFF, ell.Edwards_Affine_Arithm),
            fg.EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm),
            fg.EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm),
            fg.EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm),
            ]
        for group in curves:
            sec_grp = sfg.SecureGroup(group, group.arithm)
            b = sec_grp(group.base_pt.value)
            b = mpc.run(mpc.output(b))
            c = sec_grp(group.base_pt)
            c = mpc.run(mpc.output(c))
            self.assertEqual(b, c)
            secfld = mpc.SecFld(modulus=sec_grp.group.order)
            out = sfg.repeat_public_base_secret_output(group.base_pt, secfld(2), sec_grp)
            out = mpc.run(mpc.output(out))
            self.assertEqual(out, group.base_pt^2)
            bp4 = group.base_pt * 4
            sec_bp = sec_grp(group.base_pt)
            sec_id = sec_grp.identity
            id_out = mpc.run(mpc.output(sec_id))
            sec_bp4 = sec_bp * 4
            sec_bp4 = sec_bp4 + sec_id
            new_bp4 = mpc.run(mpc.output(sec_bp4))
            self.assertEqual(new_bp4, bp4)
            secfld = mpc.SecFld(modulus=sec_grp.group.order)
            sec_bp8 = sfg.repeat_public_base_secret_output(new_bp4, secfld(2), sec_grp)
            new_bp8 = mpc.run(mpc.output(sec_bp8))
            self.assertEqual(new_bp8, bp4 + bp4)
            sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2), sec_grp)
            new_bp8 = mpc.run(mpc.output(sec_bp8))
            self.assertEqual(new_bp8, bp4 + bp4)
            sec_bp8 = sfg.secure_repeat(new_bp4, secfld(2))
            new_bp8 = mpc.run(mpc.output(sec_bp8))
            self.assertEqual(new_bp8, bp4 + bp4)

    def test_sec_qr_group(self):
        order, modulus = find_safe_primes(64)
        group = fg.QuadraticResidue(modulus=modulus)
        sec_grp = sfg.SecureGroup(group)
        sec_g4 = sec_grp(2) * sec_grp(2)
        g4 = mpc.run(mpc.output(sec_g4))
        self.assertEqual(int(group(2)), int(group.identity * group(2)))
        self.assertEqual(g4, group(2) * group(2))
        secfld = mpc.SecFld(modulus=sec_grp.group.order)
        self.assertEqual(mpc.run(mpc.output(sfg.secure_repeat(group(2), secfld(2), sec_grp))), g4)


if __name__ == "__main__":
    unittest.main()

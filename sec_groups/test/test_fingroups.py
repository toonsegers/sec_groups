import unittest
import os
import sys

project_root = sys.path.append(os.path.abspath("../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import sec_groups.fingroups as fg
import sec_groups.ellcurves as ell
import sec_groups.pairing as pairing


class FiniteGrps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_curves(self):
        curves = [
            fg.EllipticCurve(ell.ED25519, ell.ED_AFF, ell.Edwards_Affine_Arithm),
            fg.EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm),
            fg.EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm),
            fg.EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm),
            ]
        for group in curves:
            generator = group.generator
            self.assertEqual(generator * -1, -generator)
            g4 = generator + generator + generator + generator
            self.assertEqual(g4, (generator * 4))
            self.assertEqual(generator * 0, group.identity)
            self.assertTrue(g4.on_curve())
            g2 = group.generator + group.generator
            g2 = g2 + group.identity
            g3 = g2 + group.generator
            g4 = group.generator * 4
            self.assertEqual(g4, g3 + group.generator)

            if "ED" in group.curve.name:
                eg = generator.to_extended()
                eg4 = eg + eg + eg + eg
                self.assertEqual(eg4, g4)
                self.assertEqual(eg4.to_affine(), g4.to_affine())


    def test_pairing(self):
        bn = fg.EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
        bn_twist = fg.EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)

        bn_base_pt = bn.base_pt
        bn_twist_base_pt = bn_twist.base_pt
        out = pairing.optimal_ate(bn_twist_base_pt, bn_base_pt)
        correct_result = "((21196725492326442667627602847448291017588785091108454726553939930507871042194x+43377286344385619261366633786771510260308325523330699049673437042931316243080,6519643025219345838017222895845262421933106748812739147102715685220942156402x+49167428805016148020211692433963935169227827998815966523759923767649702416286,20859153989312102139134197789193625833054352751342763983124765227328906753159x+54487609867103086000679472636440811221594782878893833473117980998935745956453),(39332187218762343173097683922586548248512461497033773500717078587710862589062x+53094021597187285527531149248961924798783924165048746910730430368152798446315,30733062774817315099333283633560206070773769397463591953591634872711994123707x+13560812839206871407210482171294630929511117297628119163695762035076749250365,57080396271372976186541856910255122494067079575964633601781325258931774013251x+60034081832369659851990276215766505463071420460830584339216728009871686767851))"
        self.assertEqual(str(out), correct_result)


    def test_symmetric_group(self):
        n = 5
        group = fg.S(n)
        self.assertEqual(group.order, 120)
        a = group.identity
        b = a @ a
        self.assertEqual(b, a)
        c = group([1, 2, 3, 4, 0])
        d = group([3, 4, 2, 1, 0])
        self.assertEqual(str(c), "(1, 2, 3, 4, 0)")
        self.assertEqual(str(d ^ n), "(3, 4, 2, 1, 0)")
        self.assertEqual(str(d ^ n + 1), "(1, 0, 2, 4, 3)")
        self.assertEqual(c != c @ a, False)

    def test_quadratic_residue_group(self):
        modulus = 101
        group = fg.QuadraticResidue(modulus)
        self.assertEqual(group.order, 50)
        a = group.identity
        b = a @ a
        self.assertEqual(b, a)
        c = group(3 ** 2)
        d = group(4 ** 2)
        self.assertEqual(group.is_cyclic, True)
        self.assertEqual(str(c.operation(d)), "43")
        self.assertEqual(str(1 / c), "45")
        self.assertEqual(str(c ** -1), "45")
        self.assertEqual(c == c * a, True)
        self.assertEqual(c != c * a, False)


if __name__ == "__main__":
    unittest.main()

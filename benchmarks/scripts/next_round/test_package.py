"""Evidence transport checks must reject changed/unmanifested/unsafe files."""
from pathlib import Path
import tempfile
import unittest
import zipfile
import numpy as np
from package import execute,digest
from verify import load_propensities


class PackageTests(unittest.TestCase):
    def fixture(self,base):
        root=Path(base)/'fp_tools_next_round_evidence';root.mkdir()
        (root/'evidence.txt').write_text('complete negative evidence\n')
        (root/'SHA256SUMS').write_text(digest(root/'evidence.txt')+'  evidence.txt\n')
        return root

    def test_verified_single_root(self):
        with tempfile.TemporaryDirectory() as base:
            root=self.fixture(base);target=Path(base)/'result.zip';execute(root,target)
            with zipfile.ZipFile(target) as z:
                self.assertEqual(len(z.namelist()),2);self.assertIsNone(z.testzip())
                self.assertTrue(all(n.startswith('fp_tools_next_round_evidence/') for n in z.namelist()))

    def test_refuse_unmanifested_or_changed(self):
        for extra in [False,True]:
            with tempfile.TemporaryDirectory() as base:
                root=self.fixture(base)
                (root/('extra.txt' if extra else 'evidence.txt')).write_text('unexpected\n')
                with self.assertRaises(ValueError):execute(root,Path(base)/'result.zip')

    def test_refuse_symlink(self):
        with tempfile.TemporaryDirectory() as base:
            root=self.fixture(base);external=Path(base)/'external';external.write_text('outside')
            (root/'linked.txt').symlink_to(external)
            with self.assertRaises(ValueError):execute(root,Path(base)/'result.zip')

    def test_shared_expectations_require_exact_site_order(self):
        with tempfile.TemporaryDirectory() as base:
            root=Path(base);(root/'arrays').mkdir()
            ids=np.array(['site_a','site_b']);expected=np.arange(402).reshape(2,201).astype(float)
            np.savez(root/'arrays/K562_rep1.propensities.npz',site_id=ids,LOG21_plus=expected,LOG21_minus=-expected)
            target=root/'arrays/K562_rep2.propensities.npz'
            np.savez(target,site_id=ids,LOG21_reference=np.array('arrays/K562_rep1.propensities.npz'))
            np.testing.assert_array_equal(load_propensities(root,'K562_rep2')['LOG21_plus'],expected)
            np.savez(target,site_id=ids[::-1],LOG21_reference=np.array('arrays/K562_rep1.propensities.npz'))
            with self.assertRaises(AssertionError):load_propensities(root,'K562_rep2')


if __name__=='__main__':unittest.main(verbosity=2)

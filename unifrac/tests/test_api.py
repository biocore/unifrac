# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, UniFrac development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

import numpy as np
import numpy.testing as npt

from unifrac import ssu


class UnifracAPITests(unittest.TestCase):
    package = 'unifrac.tests'

    def test_meta_unifrac(self):
        t1 = self.get_data_path('t1.newick')
        e1 = self.get_data_path('e1.biom')

        result = ssu(e1, t1, 'unweighted', False, 1.0, 1)

        u1_distances = np.array([[0, 10/16., 8/13.],
                                 [10/16., 0, 8/17.],
                                 [8/13., 8/17., 0]])

        npt.assert_almost_equal(u1_distances, result.data)
        self.assertEqual(tuple('ABC'), result.ids)

    def test_ssu_bad_tree(self):
        e1 = self.get_data_path('e1.biom')
        with self.assertRaisesRegex(IOError, "Tree file not found."):
            ssu(e1, 'bad-file', 'unweighted', False, 1.0, 1)

    def test_ssu_bad_table(self):
        t1 = self.get_data_path('t1.newick')
        with self.assertRaisesRegex(IOError, "Table file not found."):
            ssu('bad-file', t1, 'unweighted', False, 1.0, 1)

    def test_ssu_bad_method(self):
        t1 = self.get_data_path('t1.newick')
        e1 = self.get_data_path('e1.biom')

        with self.assertRaisesRegex(ValueError, "Unknown method."):
            ssu(e1, t1, 'unweightedfoo', False, 1.0, 1)


if __name__ == "__main__":
    unittest.main()

# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import numpy.testing as npt

from unifrac import meta


class StateUnifracTests(unittest.TestCase):
    package = 'unifrac.tests'

    def test_meta_unifrac(self):
        """meta_unifrac should give correct result on sample trees"""
        t1 = self.get_data_path('t1.newick')
        t2 = self.get_data_path('t2.newick')
        e1 = self.get_data_path('e1.biom')
        e2 = self.get_data_path('e2.biom')

        result = meta([e1, e2], [t1, t2],
                      weights=[1, 1],
                      consolidation='skipping-missing-values',
                      method='unweighted')

        u1_distances = np.array([[0, 10/16., 8/13.],
                                 [10/16., 0, 8/17.],
                                 [8/13., 8/17., 0]])
        u2_distances = np.array([[0, 11/14., 6/13.],
                                 [11/14., 0, 7/13.],
                                 [6/13., 7/13., 0]])
        exp = (u1_distances + u2_distances) / 2
        npt.assert_almost_equal(exp, result.data)
        self.assertEqual(tuple('ABC'), result.ids)

    def test_meta_unifrac_unbalanced(self):
        with self.assertRaisesRegex(ValueError, ("Number of trees and tables "
                                                 "must be the same.")):
            meta(('a', ), ('a', 'b'))

        with self.assertRaisesRegex(ValueError, ("Number of trees and tables "
                                                 "must be the same.")):
            meta(('a', 'b'), ('a', ))

    def test_meta_unifrac_unbalanced_weights(self):
        with self.assertRaisesRegex(ValueError, "Number of weights does not "
                                                "match number of trees and "
                                                "tables."):
            meta(('c', 'd'), ('a', 'b'), weights=(1, 2, 3))

    def test_meta_unifrac_missing(self):
        with self.assertRaisesRegex(ValueError, "No trees specified."):
            meta(('a', ), tuple())

        with self.assertRaisesRegex(ValueError, "No tables specified."):
            meta(tuple(), ('a', ))

    def test_meta_unifrac_no_method(self):
        with self.assertRaisesRegex(ValueError, "No method specified."):
            meta(('a', ), ('b', ))

    def test_meta_unifrac_bad_method(self):
        with self.assertRaisesRegex(ValueError, "Method \(bar\) "
                                                "unrecognized."):
            meta(('a', ), ('b', ), method='bar')

    def test_meta_unifrac_bad_consolidation(self):
        with self.assertRaisesRegex(ValueError,
                                    "Consolidation \(foo\) unrecognized."):
            meta(('a', ), ('b', ), method='unweighted', consolidation='foo')

    def test_meta_unifrac_alpha_not_generalized(self):
        with self.assertRaisesRegex(ValueError,
                                    "The alpha parameter can"):
            meta(('a', ), ('b', ), method='generalized',
                 alpha=1, consolidation='skipping_missing_matrices')


if __name__ == "__main__":
    unittest.main()

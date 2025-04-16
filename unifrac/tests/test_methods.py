# ----------------------------------------------------------------------------
# Copyright (c) 2016-, UniFrac development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import pkg_resources
import os

import h5py
import biom
import skbio
import numpy as np
import numpy.testing as npt
import pandas.testing as pdt

from unifrac import meta
from unifrac._methods import (_call_ssu, has_samples_biom_v210, has_negative,
                              has_atleast_two_samples, _call_faith_pd)


class StateUnifracTests(unittest.TestCase):
    package = 'unifrac.tests'

    def setUp(self):
        super().setUp()
        self.table1 = self.get_data_path('e1.biom')
        self.table2 = self.get_data_path('e2.biom')
        self.tree1 = self.get_data_path('t1.newick')
        self.tree2 = self.get_data_path('t2.newick')
        self.not_a_table = self.tree1
        self.not_a_tree = self.table1

    def get_data_path(self, filename):
        # adapted from qiime2.plugin.testing.TestPluginBase
        return pkg_resources.resource_filename(self.package,
                                               'data/%s' % filename)

    def test_meta_unifrac(self):
        """meta_unifrac should give correct result on sample trees"""
        result = meta([self.table1, self.table2], [self.tree1, self.tree2],
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
            meta((self.table1, ), (self.tree1, self.tree2),
                 method='unweighted')

        with self.assertRaisesRegex(ValueError, ("Number of trees and tables "
                                                 "must be the same.")):
            meta((self.table1, self.table2), (self.tree1, ),
                 method='unweighted')

    def test_meta_unifrac_unbalanced_weights(self):
        with self.assertRaisesRegex(ValueError, "Number of weights does not "
                                                "match number of trees and "
                                                "tables."):
            meta((self.table1, self.table2), (self.tree1, self.tree2),
                 weights=(1, 2, 3), )

    def test_meta_unifrac_missing(self):
        with self.assertRaisesRegex(ValueError, "No trees specified."):
            meta((self.table1, ), tuple(), method='unweighted')

        with self.assertRaisesRegex(ValueError, "No tables specified."):
            meta(tuple(), (self.tree1, ), method='unweighted')

    def test_meta_validation(self):
        with self.assertRaisesRegex(ValueError,
                                    "Table does not appear to be a "
                                    "BIOM-Format v2.1"):
            meta((self.table1, self.not_a_table), (self.tree1, self.tree2),
                 method='unweighted')

        with self.assertRaisesRegex(ValueError,
                                    "The phylogeny does not "
                                    "appear to be newick"):
            meta((self.table1, self.table2), (self.tree1, self.not_a_tree),
                 method='unweighted')

    def test_meta_unifrac_no_method(self):
        with self.assertRaisesRegex(ValueError, "No method specified."):
            meta((self.table1, ), (self.tree1, ))

    def test_meta_unifrac_bad_method(self):
        with self.assertRaisesRegex(ValueError, r"Method \(bar\) "
                                                "unrecognized."):
            meta((self.table1, ), (self.tree1, ), method='bar')

    def test_meta_unifrac_bad_consolidation(self):
        with self.assertRaisesRegex(ValueError,
                                    r"Consolidation \(foo\) unrecognized."):
            meta((self.table1, ), (self.tree1, ), method='unweighted',
                 consolidation='foo')

    def test_meta_unifrac_alpha_not_generalized(self):
        with self.assertRaisesRegex(ValueError,
                                    "The alpha parameter can"):
            meta((self.table1, ), (self.tree1, ), method='unweighted',
                 alpha=1, consolidation='skipping_missing_matrices')

    def test_has_samples_biom_v210(self):
        fp = self.get_data_path('crawford.biom')
        self.assertTrue(has_samples_biom_v210(fp))

        tmpfile = '/tmp/fake.biom'
        empty = biom.Table([], [], [])
        try:
            with h5py.File(tmpfile, 'w') as fp:
                empty.to_hdf5(fp, 'asd')
            fp = self.get_data_path(tmpfile)
            self.assertFalse(has_samples_biom_v210(tmpfile))
        finally:
            os.unlink(tmpfile)

    def test_has_negative(self):
        self.assertFalse(has_negative(self.get_data_path('crawford.biom')))
        tab = biom.load_table(self.get_data_path('crawford.biom'))
        tab._data[1] *= -1
        tmpfile = '/tmp/fake.biom'
        try:
            with h5py.File(tmpfile, 'w') as fp:
                tab.to_hdf5(fp, 'asd')
            fp = self.get_data_path(tmpfile)
            self.assertTrue(has_negative(tmpfile))
        finally:
            os.unlink(tmpfile)

    def test_has_atleast_two_samples(self):
        path = self.get_data_path('crawford.biom')
        self.assertTrue(has_atleast_two_samples(path))

        tab = biom.Table([], [], [])
        tmpfile = '/tmp/fake.biom'
        try:
            with h5py.File(tmpfile, 'w') as fp:
                tab.to_hdf5(fp, 'asd')
            fp = self.get_data_path(tmpfile)
            self.assertFalse(has_atleast_two_samples(tmpfile))
        finally:
            os.unlink(tmpfile)

    def test_call_ssu_empty_biom(self):
        empty = biom.Table([], [], [])
        tre = skbio.TreeNode()
        with self.assertRaisesRegex(ValueError, "contain any samples"):
            _call_ssu(empty, tre)

    def test_call_faith_pd(self):
        obs_1 = _call_faith_pd(self.get_data_path('crawford.biom'),
                               self.get_data_path('crawford.tre'))

        tab = biom.load_table(self.get_data_path('crawford.biom'))
        tre = skbio.TreeNode.read(self.get_data_path('crawford.tre'))

        obs_2 = _call_faith_pd(tab, tre)

        pdt.assert_series_equal(obs_1, obs_2)
        self.assertEqual(len(tab.ids()), len(obs_1))


if __name__ == "__main__":
    unittest.main()

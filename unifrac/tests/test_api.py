# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, UniFrac development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import os
from io import StringIO
from tempfile import gettempdir
import pkg_resources

import numpy as np
import numpy.testing as npt
from biom import Table, load_table
from biom.util import biom_open
from skbio import TreeNode
import skbio.diversity

from unifrac import ssu, faith_pd


class UnifracAPITests(unittest.TestCase):
    package = 'unifrac.tests'

    def get_data_path(self, filename):
        # adapted from qiime2.plugin.testing.TestPluginBase
        return pkg_resources.resource_filename(self.package,
                                               'data/%s' % filename)

    def test_unweighted_root_eval_issue_46(self):
        tree = self.get_data_path('crawford.tre')
        table = self.get_data_path('crawford.biom')

        table_inmem = load_table(table)
        tree_inmem = skbio.TreeNode.read(tree)

        ids = table_inmem.ids()
        otu_ids = table_inmem.ids(axis='observation')
        cnts = table_inmem.matrix_data.astype(int).toarray().T
        exp = skbio.diversity.beta_diversity('unweighted_unifrac', cnts,
                                             ids=ids, otu_ids=otu_ids,
                                             tree=tree_inmem)
        obs = ssu(table, tree, 'unweighted', False, 1.0, False, 1)
        npt.assert_almost_equal(obs.data, exp.data)

    def test_meta_unifrac(self):
        t1 = self.get_data_path('t1.newick')
        e1 = self.get_data_path('e1.biom')

        result = ssu(e1, t1, 'unweighted', False, 1.0, False, 1)

        u1_distances = np.array([[0, 10 / 16., 8 / 13.],
                                 [10 / 16., 0, 8 / 17.],
                                 [8 / 13., 8 / 17., 0]])

        npt.assert_almost_equal(u1_distances, result.data)
        self.assertEqual(tuple('ABC'), result.ids)

    def test_ssu_bad_tree(self):
        e1 = self.get_data_path('e1.biom')
        with self.assertRaisesRegex(IOError, "Tree file not found."):
            ssu(e1, 'bad-file', 'unweighted', False, 1.0, False, 1)

    def test_ssu_bad_table(self):
        t1 = self.get_data_path('t1.newick')
        with self.assertRaisesRegex(IOError, "Table file not found."):
            ssu('bad-file', t1, 'unweighted', False, 1.0, False, 1)

    def test_ssu_bad_method(self):
        t1 = self.get_data_path('t1.newick')
        e1 = self.get_data_path('e1.biom')

        with self.assertRaisesRegex(ValueError, "Unknown method."):
            ssu(e1, t1, 'unweightedfoo', False, 1.0, False, 1)


class EdgeCasesTests(unittest.TestCase):
    # These tests were mostly ported from skbio's
    # skbio/diversity/beta/tests/test_unifrac.py at SHA-256 ea901b3b6b0b
    # note that not all tests were kept since the APIs are different.
    #
    # The test cases below only exercise unweighted, weighted and weighted
    # normalized UniFrac. The C++ test suite verifies (against reference
    # implementations) the variance adjusted and generalized variants of the
    # algorithm.

    package = 'unifrac.tests'

    def _work(self, u_counts, v_counts, otu_ids, tree, method):
        data = np.array([u_counts, v_counts]).T

        bt = Table(data, otu_ids, ['u', 'v'])

        ta = os.path.join(gettempdir(), 'table.biom')
        tr = os.path.join(gettempdir(), 'tree.biom')

        self.files_to_delete.append(ta)
        self.files_to_delete.append(tr)

        with biom_open(ta, 'w') as fhdf5:
            bt.to_hdf5(fhdf5, 'Table for unit testing')
        tree.write(tr)

        # return value is a distance matrix, get the distance from u->v
        return ssu(ta, tr, method, False, 1.0, False, 1)['u', 'v']

    def weighted_unifrac(self, u_counts, v_counts, otu_ids, tree,
                         normalized=False):
        if normalized:
            method = 'weighted_normalized'
        else:
            method = 'weighted_unnormalized'
        return self._work(u_counts, v_counts, otu_ids, tree, method)

    def unweighted_unifrac(self, u_counts, v_counts, otu_ids, tree,
                           normalized=False):
        return self._work(u_counts, v_counts, otu_ids, tree, 'unweighted')

    def setUp(self):
        self.b1 = np.array(
            [[1, 3, 0, 1, 0],
             [0, 2, 0, 4, 4],
             [0, 0, 6, 2, 1],
             [0, 0, 1, 1, 1],
             [5, 3, 5, 0, 0],
             [0, 0, 0, 3, 5]])
        self.sids1 = list('ABCDEF')
        self.oids1 = ['OTU%d' % i for i in range(1, 6)]
        self.t1 = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        self.t1_w_extra_tips = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,(OTU5:0.25,(OTU6:0.5,OTU7:0.5):0.5):0.5):1.25):0.0'
                     ')root;'))

        self.t2 = TreeNode.read(
            StringIO('((OTU1:0.1, OTU2:0.2):0.3, (OTU3:0.5, OTU4:0.7):1.1)'
                     'root;'))
        self.oids2 = ['OTU%d' % i for i in range(1, 5)]

        self.files_to_delete = []

    def tearDown(self):

        for f in self.files_to_delete:
            try:
                os.remove(f)
            except OSError:
                pass

    def test_ssu_table_not_subset_tree(self):
        tree = TreeNode.read(StringIO('((OTU1:0.5,OTU3:1.0):1.0)root;'))
        expected_message = "The table does not appear to be completely "\
                           "represented by the phylogeny."

        with self.assertRaisesRegex(ValueError, expected_message):
            self.unweighted_unifrac(self.b1[0], self.b1[1], self.oids1, tree)

    def test_unweighted_otus_out_of_order(self):
        # UniFrac API does not assert the observations are in tip order of the
        # input tree
        shuffled_ids = self.oids1[:]
        shuffled_b1 = self.b1.copy()

        shuffled_ids[0], shuffled_ids[-1] = shuffled_ids[-1], shuffled_ids[0]
        shuffled_b1[:, [0, -1]] = shuffled_b1[:, [-1, 0]]

        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = self.unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = self.unweighted_unifrac(
                    shuffled_b1[i], shuffled_b1[j], shuffled_ids, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_otus_out_of_order(self):
        # UniFrac API does not assert the observations are in tip order of the
        # input tree
        shuffled_ids = self.oids1[:]
        shuffled_b1 = self.b1.copy()

        shuffled_ids[0], shuffled_ids[-1] = shuffled_ids[-1], shuffled_ids[0]
        shuffled_b1[:, [0, -1]] = shuffled_b1[:, [-1, 0]]

        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = self.weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = self.weighted_unifrac(
                    shuffled_b1[i], shuffled_b1[j], shuffled_ids, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_unweighted_extra_tips(self):
        # UniFrac values are the same despite unobserved tips in the tree
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = self.unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1_w_extra_tips)
                expected = self.unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_extra_tips(self):
        # UniFrac values are the same despite unobserved tips in the tree
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = self.weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1_w_extra_tips)
                expected = self.weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_unweighted_minimal_trees(self):
        # two tips
        tree = TreeNode.read(StringIO('(OTU1:0.25, OTU2:0.25)root;'))
        actual = self.unweighted_unifrac([1, 0], [0, 0], ['OTU1', 'OTU2'],
                                         tree)
        expected = 1.0
        self.assertEqual(actual, expected)

    def test_unweighted_root_not_observed(self):
        # expected values computed with QIIME 1.9.1 and by hand
        # root node not observed, but branch between (OTU1, OTU2) and root
        # is considered shared
        actual = self.unweighted_unifrac([1, 1, 0, 0], [1, 0, 0, 0],
                                         self.oids2, self.t2)
        # for clarity of what I'm testing, compute expected as it would
        # based on the branch lengths. the values that compose shared was
        # a point of confusion for me here, so leaving these in for
        # future reference
        expected = 0.2 / (0.1 + 0.2 + 0.3)  # 0.3333333333
        self.assertAlmostEqual(actual, expected)

        # root node not observed, but branch between (OTU3, OTU4) and root
        # is considered shared
        actual = self.unweighted_unifrac([0, 0, 1, 1], [0, 0, 1, 0],
                                         self.oids2, self.t2)
        # for clarity of what I'm testing, compute expected as it would
        # based on the branch lengths. the values that compose shared was
        # a point of confusion for me here, so leaving these in for
        # future reference
        expected = 0.7 / (1.1 + 0.5 + 0.7)  # 0.3043478261
        self.assertAlmostEqual(actual, expected)

    def test_weighted_root_not_observed(self):
        # expected values computed by hand, these disagree with QIIME 1.9.1
        # root node not observed, but branch between (OTU1, OTU2) and root
        # is considered shared
        actual = self.weighted_unifrac([1, 0, 0, 0], [1, 1, 0, 0],
                                       self.oids2, self.t2)
        expected = 0.15
        self.assertAlmostEqual(actual, expected)

        # root node not observed, but branch between (OTU3, OTU4) and root
        # is considered shared
        actual = self.weighted_unifrac([0, 0, 1, 1], [0, 0, 1, 0],
                                       self.oids2, self.t2)
        expected = 0.6
        self.assertAlmostEqual(actual, expected)

    def test_weighted_normalized_root_not_observed(self):
        # expected values computed by hand, these disagree with QIIME 1.9.1
        # root node not observed, but branch between (OTU1, OTU2) and root
        # is considered shared
        actual = self.weighted_unifrac([1, 0, 0, 0], [1, 1, 0, 0],
                                       self.oids2, self.t2, normalized=True)
        expected = 0.1764705882
        self.assertAlmostEqual(actual, expected)

        # root node not observed, but branch between (OTU3, OTU4) and root
        # is considered shared
        actual = self.weighted_unifrac([0, 0, 1, 1], [0, 0, 1, 0],
                                       self.oids2, self.t2, normalized=True)
        expected = 0.1818181818
        self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac_identity(self):
        for i in range(len(self.b1)):
            actual = self.unweighted_unifrac(
                self.b1[i], self.b1[i], self.oids1, self.t1)
            expected = 0.0
            self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac_symmetry(self):
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = self.unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = self.unweighted_unifrac(
                    self.b1[j], self.b1[i], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac_non_overlapping(self):
        # these communities only share the root node
        actual = self.unweighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 1, 1], self.oids1, self.t1)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # unweighted unifrac implementation
        # sample A versus all
        actual = self.unweighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1)
        expected = 0.238095238095
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[0], self.b1[2], self.oids1, self.t1)
        expected = 0.52
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[0], self.b1[3], self.oids1, self.t1)
        expected = 0.52
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[0], self.b1[4], self.oids1, self.t1)
        expected = 0.545454545455
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[0], self.b1[5], self.oids1, self.t1)
        expected = 0.619047619048
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = self.unweighted_unifrac(
            self.b1[1], self.b1[2], self.oids1, self.t1)
        expected = 0.347826086957
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[1], self.b1[3], self.oids1, self.t1)
        expected = 0.347826086957
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[1], self.b1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[1], self.b1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = self.unweighted_unifrac(
            self.b1[2], self.b1[3], self.oids1, self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[2], self.b1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[2], self.b1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = self.unweighted_unifrac(
            self.b1[3], self.b1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = self.unweighted_unifrac(
            self.b1[3], self.b1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = self.unweighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_identity(self):
        for i in range(len(self.b1)):
            actual = self.weighted_unifrac(
                self.b1[i], self.b1[i], self.oids1, self.t1)
            expected = 0.0
            self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_symmetry(self):
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = self.weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = self.weighted_unifrac(
                    self.b1[j], self.b1[i], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_non_overlapping(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # weighted unifrac implementation
        # these communities only share the root node
        actual = self.weighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 4.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # weighted unifrac implementation
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1)
        expected = 2.4
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[2], self.oids1, self.t1)
        expected = 1.86666666667
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[3], self.oids1, self.t1)
        expected = 2.53333333333
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[4], self.oids1, self.t1)
        expected = 1.35384615385
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[5], self.oids1, self.t1)
        expected = 3.2
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = self.weighted_unifrac(
            self.b1[1], self.b1[2], self.oids1, self.t1)
        expected = 2.26666666667
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[1], self.b1[3], self.oids1, self.t1)
        expected = 0.933333333333
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[1], self.b1[4], self.oids1, self.t1)
        expected = 3.2
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[1], self.b1[5], self.oids1, self.t1)
        expected = 0.8375
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = self.weighted_unifrac(
            self.b1[2], self.b1[3], self.oids1, self.t1)
        expected = 1.33333333333
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[2], self.b1[4], self.oids1, self.t1)
        expected = 1.89743589744
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[2], self.b1[5], self.oids1, self.t1)
        expected = 2.66666666667
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = self.weighted_unifrac(
            self.b1[3], self.b1[4], self.oids1, self.t1)
        expected = 2.66666666667
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[3], self.b1[5], self.oids1, self.t1)
        expected = 1.33333333333
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = self.weighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 4.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_identity_normalized(self):
        for i in range(len(self.b1)):
            actual = self.weighted_unifrac(
                self.b1[i], self.b1[i], self.oids1, self.t1, normalized=True)
            expected = 0.0
            self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_symmetry_normalized(self):
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = self.weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1,
                    normalized=True)
                expected = self.weighted_unifrac(
                    self.b1[j], self.b1[i], self.oids1, self.t1,
                    normalized=True)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_non_overlapping_normalized(self):
        # these communities only share the root node
        actual = self.weighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 1, 1], self.oids1, self.t1,
            normalized=True)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_normalized(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # weighted unifrac implementation
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1, normalized=True)
        expected = 0.6
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[2], self.oids1, self.t1, normalized=True)
        expected = 0.466666666667
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[3], self.oids1, self.t1, normalized=True)
        expected = 0.633333333333
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[4], self.oids1, self.t1, normalized=True)
        expected = 0.338461538462
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[0], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 0.8
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = self.weighted_unifrac(
            self.b1[1], self.b1[2], self.oids1, self.t1, normalized=True)
        expected = 0.566666666667
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[1], self.b1[3], self.oids1, self.t1, normalized=True)
        expected = 0.233333333333
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[1], self.b1[4], self.oids1, self.t1, normalized=True)
        expected = 0.8
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[1], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 0.209375
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = self.weighted_unifrac(
            self.b1[2], self.b1[3], self.oids1, self.t1, normalized=True)
        expected = 0.333333333333
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[2], self.b1[4], self.oids1, self.t1, normalized=True)
        expected = 0.474358974359
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[2], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 0.666666666667
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = self.weighted_unifrac(
            self.b1[3], self.b1[4], self.oids1, self.t1, normalized=True)
        expected = 0.666666666667
        self.assertAlmostEqual(actual, expected)
        actual = self.weighted_unifrac(
            self.b1[3], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 0.333333333333
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = self.weighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)


class FaithPDEdgeCasesTests(unittest.TestCase):
    # These tests were mostly ported from skbio's
    # skbio/diversity/alpha/tests/test_fatih_pd.py at SHA-256 a8c086b
    # note that not all tests were kept since the APIs are different.

    package = 'unifrac.tests'

    def write_table_tree(self, u_counts, otu_ids, sample_ids, tree):
        data = np.array([u_counts]).T

        bt = Table(data, otu_ids, sample_ids)

        ta = os.path.join(gettempdir(), 'table.biom')
        tr = os.path.join(gettempdir(), 'tree.biom')

        self.files_to_delete.append(ta)
        self.files_to_delete.append(tr)

        with biom_open(ta, 'w') as fhdf5:
            bt.to_hdf5(fhdf5, 'Table for unit testing')
        tree.write(tr)

        return ta, tr

    def faith_pd_work(self, u_counts, otu_ids, sample_ids, tree):
        ta, tr = self.write_table_tree(u_counts, otu_ids, sample_ids, tree)

        return faith_pd(ta, tr)

    def setUp(self):
        self.counts = np.array([0, 1, 1, 4, 2, 5, 2, 4, 1, 2])
        self.b1 = np.array([[1, 3, 0, 1, 0],
                            [0, 2, 0, 4, 4],
                            [0, 0, 6, 2, 1],
                            [0, 0, 1, 1, 1]])
        self.sids1 = list('ABCD')
        self.oids1 = ['OTU%d' % i for i in range(1, 6)]
        self.t1 = TreeNode.read(StringIO(
            '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):'
            '0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;'))
        self.t1_w_extra_tips = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,(OTU5:0.25,(OTU6:0.5,OTU7:0.5):0.5):0.5):1.25):0.0'
                     ')root;'))

        self.files_to_delete = []

    def tearDown(self):

        for f in self.files_to_delete:
            try:
                os.remove(f)
            except OSError:
                pass

    def test_faith_pd_zero_branches_omitted(self):
        # also deleted branch length fo

        t2 = TreeNode.read(StringIO(
            '((OTU1:0.5,OTU2:0.5),(OTU3:1.0,(OTU4:0.5,'
            'OTU5:0.75):1.0):1.0)root;'
        ))

        actual = self.faith_pd_work([1, 1, 0, 0, 0], self.oids1, ['foo'], t2)
        expected = 1.0
        self.assertAlmostEqual(actual[0], expected)

    def test_faith_pd_none_observed(self):
        actual = self.faith_pd_work([0, 0, 0, 0, 0], self.oids1, ['foo'],
                                    self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual.values, expected)

    def test_faith_pd_biom_table_empty(self):
        table, tree = self.write_table_tree([], [], [],
                                            self.t1)

        self.assertRaises(ValueError, faith_pd, table, tree)

    def test_faith_pd_table_not_subset_tree(self):
        tree = TreeNode.read(StringIO('((OTU1:0.5,OTU3:1.0):1.0)root;'))
        table_ids = ['OTU1', 'OTU2']
        table, tree = self.write_table_tree([1, 0], table_ids, ['foo'],
                                            tree)
        expected_message = "The table does not appear to be completely "\
                           "represented by the phylogeny."
        with self.assertRaisesRegex(ValueError, expected_message):
            faith_pd(table, tree)

    def test_faith_pd_all_observed(self):
        actual = self.faith_pd_work([1, 1, 1, 1, 1], self.oids1, ['foo'],
                                    self.t1)
        expected = sum(n.length for n in self.t1.traverse()
                       if n.length is not None)
        self.assertAlmostEqual(actual.values, expected)

        actual = self.faith_pd_work([1, 2, 3, 4, 5], self.oids1, ['foo'],
                                    self.t1)
        expected = sum(n.length for n in self.t1.traverse()
                       if n.length is not None)
        self.assertAlmostEqual(actual.values, expected)

    def test_faith_pd(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation unifrac's initial
        # phylogenetic diversity implementation
        actual = self.faith_pd_work(self.b1[0], self.oids1, [self.sids1[0]],
                                    self.t1)
        expected = 4.5
        self.assertAlmostEqual(actual.values, expected)
        actual = self.faith_pd_work(self.b1[1], self.oids1, [self.sids1[1]],
                                    self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual.values, expected)
        actual = self.faith_pd_work(self.b1[2], self.oids1, [self.sids1[2]],
                                    self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual.values, expected)
        actual = self.faith_pd_work(self.b1[3], self.oids1, [self.sids1[3]],
                                    self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual.values, expected)

    def test_faith_pd_extra_tips(self):
        # results are the same despite presences of unobserved tips in tree
        actual = self.faith_pd_work(self.b1[0], self.oids1, [self.sids1[0]],
                                    self.t1_w_extra_tips)
        expected = self.faith_pd_work(self.b1[0], self.oids1, [self.sids1[0]],
                                      self.t1)
        self.assertAlmostEqual(actual.values, expected.values)
        actual = self.faith_pd_work(self.b1[1], self.oids1, [self.sids1[1]],
                                    self.t1_w_extra_tips)
        expected = self.faith_pd_work(self.b1[1], self.oids1, [self.sids1[1]],
                                      self.t1)
        self.assertAlmostEqual(actual.values, expected.values)
        actual = self.faith_pd_work(self.b1[2], self.oids1, [self.sids1[2]],
                                    self.t1_w_extra_tips)
        expected = self.faith_pd_work(self.b1[2], self.oids1, [self.sids1[2]],
                                      self.t1)
        self.assertAlmostEqual(actual.values, expected.values)
        actual = self.faith_pd_work(self.b1[3], self.oids1, [self.sids1[3]],
                                    self.t1_w_extra_tips)
        expected = self.faith_pd_work(self.b1[3], self.oids1, [self.sids1[3]],
                                      self.t1)
        self.assertAlmostEqual(actual.values, expected.values)

    def test_faith_pd_minimal(self):
        # two tips
        tree = TreeNode.read(StringIO('(OTU1:0.25, OTU2:0.25)root;'))
        actual = self.faith_pd_work([1, 0], ['OTU1', 'OTU2'], ['foo'], tree)
        expected = 0.25
        self.assertEqual(actual.values, expected)

    def test_faith_pd_series_name(self):
        tree = TreeNode.read(StringIO('(OTU1:0.25, OTU2:0.25)root;'))
        actual = self.faith_pd_work([1, 0], ['OTU1', 'OTU2'], ['foo'], tree)
        self.assertEqual("faith_pd", actual.name)

    def test_faith_pd_root_not_observed(self):
        # expected values computed by hand
        tree = TreeNode.read(
            StringIO('((OTU1:0.1, OTU2:0.2):0.3, (OTU3:0.5, OTU4:0.7):1.1)'
                     'root;'))
        otu_ids = ['OTU%d' % i for i in range(1, 5)]
        # root node not observed, but branch between (OTU1, OTU2) and root
        # is considered observed
        actual = self.faith_pd_work([1, 1, 0, 0], otu_ids, ['foo'], tree)
        expected = 0.6
        self.assertAlmostEqual(actual[0], expected)

        # root node not observed, but branch between (OTU3, OTU4) and root
        # is considered observed
        actual = self.faith_pd_work([0, 0, 1, 1], otu_ids, ['foo'], tree)
        expected = 2.3
        self.assertAlmostEqual(actual[0], expected)

    def test_faith_pd_invalid_input(self):
        # tests are based of skbio tests, checking for duplicate ids,
        # negative counts are not included but should be incorporated

        # tree has duplicated tip ids
        tree = TreeNode.read(
            StringIO('((OTU1:0.1, OTU2:0.2):0.3, (OTU3:0.5, OTU4:0.7):1.1)'
                     'root;'))
        otu_ids = ['OTU%d' % i for i in range(1, 5)]
        u_counts = [1, 1, 0, 0]

        data = np.array([u_counts]).T

        bt = Table(data, otu_ids, ['u'])

        ta = os.path.join(gettempdir(), 'table.biom')
        tr = os.path.join(gettempdir(), 'tree.biom')

        self.files_to_delete.append(ta)
        self.files_to_delete.append(tr)

        with biom_open(ta, 'w') as fhdf5:
            bt.to_hdf5(fhdf5, 'Table for unit testing')
        tree.write(tr)

        self.assertRaises(IOError, faith_pd, 'dne.biom', tr)
        self.assertRaises(IOError, faith_pd, ta, 'dne.tre')


if __name__ == "__main__":
    unittest.main()

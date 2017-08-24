# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import skbio
from q2_types.feature_table import BIOMV210Format
from q2_types.tree import NewickFormat
from q2_state_unifrac._meta import CONSOLIDATIONS

from q2_state_unifrac._api import ssu


def unweighted(table: BIOMV210Format,
               phylogeny: NewickFormat,
               threads: int=1,
               variance_adjusted: bool=False)-> skbio.DistanceMatrix:
    return ssu(str(table), str(tree), output_fp, 'unweighted',
               variance_adjusted, 1.0, threads)


def weighted_normalized(table: BIOMV210Format,
                        phylogeny: NewickFormat,
                        threads: int=1,
                        variance_adjusted: bool=False)-> skbio.DistanceMatrix:
    return ssu(str(table), str(tree), output_fp, 'weighted_normalized',
               variance_adjusted, 1.0, threads)


def weighted_unnormalized(table: BIOMV210Format,
                          phylogeny: NewickFormat,
                          threads: int=1,
                          variance_adjusted: bool=False) -> skbio.DistanceMatrix:  # noqa
    return ssu(str(table), str(tree), output_fp, 'weighted_unnormalized',
               variance_adjusted, 1.0, threads)


def generalized(table: BIOMV210Format,
                phylogeny: NewickFormat,
                threads: int=1,
                alpha: float=1.0,
                variance_adjusted: bool=False)-> skbio.DistanceMatrix:
    if alpha == 1.0:
        warn("alpha of 1.0 is weighted-normalized UniFrac. "
             "Weighted-normalized is being used instead as it is more "
             "optimized.",
             Warning)
        return weighted_normalized(table, phylogeny, threads,
                                   variance_adjusted)
    else:
        return ssu(str(table), str(tree), output_fp, 'generalized',
                   variance_adjusted, alpha, threads)

METHODS = {'unweighted': unweighted,
           'weighted_normalized': weighted_normalized,
           'weighted_unnormalized': weighted_unnormalized,
           'generalized': generalized}


def meta(tables: tuple, phylogenies: tuple, weights: tuple=None,
         consolidation: str=None, method: str=None,
         threads: int=1, variance_adjusted: bool=False,
         alpha: float=None) -> skbio.DistanceMatrix:
    if not len(tables):
        raise ValueError("No tables specified.")

    if not len(phylogenies):
        raise ValueError("No trees specified.")

    if len(tables) != len(phylogenies):
        raise ValueError("Number of trees and tables must be the same.")

    if weights is None:
        weights = tuple(1 for _ in phylogenies)
    else:
        if len(weights) != len(phylogenies):
            raise ValueError("Number of weights does not match number of "
                             "trees and tables.")

    if method is None:
        raise ValueError("No method specified.")
    method_ = METHODS.get(method.replace('-', '_'))
    if method_ is None:
        raise ValueError("Method (%s) unrecognized. Available methods are: %s"
                         % (method, ', '.join(METHODS.keys())))

    if consolidation is None:
        raise ValueError("No consolidation specified.")
    consolidation_ = CONSOLIDATIONS.get(consolidation.replace('-', '_'))
    if consolidation_ is None:
        raise ValueError("Consolidation (%s) unrecognized. Available "
                         "consolidations are: %s"
                         % (consolidation, ', '.join(CONSOLIDATIONS.keys())))

    if alpha is not None and method is not generalized:
        raise ValueError("The alpha parameter can only be set when the method "
                         "is set as 'generalized', the selected method is "
                         "'%s'." % method)

    kwargs = {'threads': threads, 'variance_adjusted': variance_adjusted}
    if alpha is not None:
        kwargs['alpha'] = alpha

    weights = np.array(weights, float)/sum(weights)
    dms = [method_(table, tree, **kwargs) for table, tree in zip(tables,
                                                                 phylogenies)]
    all_ids = sorted(reduce(or_, [set(dm.ids) for dm in dms]))
    dm = consolidation_(dms, [dm.ids for dm in dms], weights, all_ids)
    return skbio.DistanceMatrix(dm, ids=all_ids)

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
    return ssu(str(table), str(tree), output_fp, 'generalized',
               variance_adjusted, alpha, threads)

# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
import tempfile

import skbio
from q2_types.feature_table import BIOMV210Format
from q2_types.tree import NewickFormat

from pkg_resources import resource_exists, Requirement, resource_filename

ARGS = (Requirement.parse('q2_state_unifrac'), 'q2_state_unifrac/ssu')


def _sanity():
    if not resource_exists(*ARGS):
        raise ValueError("ssu could not be located!")


def _run(table_fp, tree_fp, output_fp, threads, method,
         variance_adjusted=False, alpha=None):
    cmd = [resource_filename(*ARGS),
           '-i', table_fp,
           '-t', tree_fp,
           '-o', output_fp,
           '-n', threads,
           '-m', method]

    if variance_adjusted:
        cmd.append('--vaw')

    if alpha is not None:
        cmd.append('-a')
        cmd.append(str(alpha))

    subprocess.run(cmd, check=True)


def unweighted(table: BIOMV210Format,
               phylogeny: NewickFormat,
               threads: int=1,
               variance_adjusted: bool=False)-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        output_fp = os.path.join(tmp, 'foo.dm')
        _run(str(table), str(phylogeny), output_fp, str(threads),
             'unweighted', variance_adjusted)
        return skbio.DistanceMatrix.read(output_fp)


def weighted_normalized(table: BIOMV210Format,
                        phylogeny: NewickFormat,
                        threads: int=1,
                        variance_adjusted: bool=False)-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        output_fp = os.path.join(tmp, 'foo.dm')
        _run(str(table), str(phylogeny), output_fp, str(threads),
             'weighted_normalized', variance_adjusted)
        return skbio.DistanceMatrix.read(output_fp)


def weighted_unnormalized(table: BIOMV210Format,
                          phylogeny: NewickFormat,
                          threads: int=1,
                          variance_adjusted: bool=False) -> skbio.DistanceMatrix:  # noqa
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        output_fp = os.path.join(tmp, 'foo.dm')
        _run(str(table), str(phylogeny), output_fp, str(threads),
             'weighted_unnormalized', variance_adjusted)
        return skbio.DistanceMatrix.read(output_fp)


def generalized(table: BIOMV210Format,
                phylogeny: NewickFormat,
                threads: int=1,
                alpha: float=1.0,
                variance_adjusted: bool=False)-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        output_fp = os.path.join(tmp, 'foo.dm')
        _run(str(table), str(phylogeny), output_fp, str(threads),
             'generalized', variance_adjusted, alpha)
        return skbio.DistanceMatrix.read(output_fp)

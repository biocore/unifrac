# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
import shutil
import tempfile

import skbio
from q2_types.feature_table import BIOMV210Format
from q2_types.tree import NewickFormat


def _sanity():
    if shutil.which('ssu') is None:
        raise ValueError("ssu does not appear in $PATH")


def _run(table_fp, tree_fp, output_fp, threads, method):
    cmd = ['ssu',
           '-i', table_fp,
           '-t', tree_fp,
           '-o', output_fp,
           '-n', threads,
           '-m', method]

    subprocess.run(cmd, check=True)


def unweighted(table: BIOMV210Format,
               phylogeny: NewickFormat,
               threads: int=1)-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        output_fp = os.path.join(tmp, 'foo.dm')
        _run(str(table), str(phylogeny), output_fp, str(threads),
             'unweighted')
        return skbio.DistanceMatrix.read(output_fp)


def weighted_normalized(table: BIOMV210Format,
                        phylogeny: NewickFormat,
                        threads: int=1)-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        output_fp = os.path.join(tmp, 'foo.dm')
        _run(str(table), str(phylogeny), output_fp, str(threads),
             'weighted_normalized')
        return skbio.DistanceMatrix.read(output_fp)


def weighted_unnormalized(table: BIOMV210Format,
                          phylogeny: NewickFormat,
                          threads: int=1)-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        output_fp = os.path.join(tmp, 'foo.dm')
        _run(str(table), str(phylogeny), output_fp, str(threads),
             'weighted_unnormalized')
        return skbio.DistanceMatrix.read(output_fp)

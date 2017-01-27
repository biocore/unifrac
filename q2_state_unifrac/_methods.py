import subprocess
import shutil
import tempfile

import h5py
import biom
import skbio


def _sanity():
    if shutil('su') is None:
        raise ValueError("su does not appear in $PATH")


def _stage(tmp, table, tree):
    table_fp = os.path.join(tmp, 'foo.biom')
    tree_fp = os.path.join(tmp, 'foo.tree')
    output_fp = os.path.join(tmp, 'foo.dm')

    with h5py.File(table_fp, 'w') as fp:
        table.to_hdf5(fp, 'q2-state-unifrac')
    with open(tree_fp, 'w') as fp:
        fp.write(str(tree))

    return (table_fp, tree_fp, output_fp)


def _run(table_fp, tree_fp, output_fp, method):
    cmd = ['su',
           '-i %s' % table_fp,
           '-t %s' % tree_fp,
           '-o %s' % output_fp,
           '-m %s' % method]

    subprocess.run(cmd, check=True)


def unweighted(table: biom.Table,
               phylogeny: skbio.TreeNode )-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        table_fp, tree_fp, output_fp = stage(tmp, table, phylogeny)
        _run(table_fp, tree_fp, output_fp, 'unweighted')
        return skbio.DistanceMatrix.read(output_fp)


def weighted_normalized(table: biom.Table,
                        phylogeny: skbio.TreeNode )-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        table_fp, tree_fp, output_fp = stage(tmp, table, phylogeny)
        _run(table_fp, tree_fp, output_fp, 'weighted_normalized')
        return skbio.DistanceMatrix.read(output_fp)


def weighted_unnormalized_(table: biom.Table,
                           phylogeny: skbio.TreeNode )-> skbio.DistanceMatrix:
    _sanity()

    with tempfile.TemporaryDirectory() as tmp:
        table_fp, tree_fp, output_fp = stage(tmp, table, phylogeny)
        _run(table_fp, tree_fp, output_fp, 'weighted_unnormalized')
        return skbio.DistanceMatrix.read(output_fp)

import biom
import skbio
import sys
import h5py

tree = skbio.TreeNode.read(sys.argv[1])
table = biom.load_table(sys.argv[2])
depth = int(sys.argv[3])

table = table.subsample(1000)
table = table.subsample(depth, by_id=True)
tree = tree.shear(set(table.ids(axis='observation')))

with h5py.File(sys.argv[4] + '/bench_%d.biom' % depth, 'w') as fp:
    table.to_hdf5(fp, 'asd')

tree.write(sys.argv[4] + '/bench_%d.tre' % depth)

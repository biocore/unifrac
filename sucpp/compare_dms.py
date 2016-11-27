import skbio
import sys
import numpy as np

dm1 = skbio.DistanceMatrix.read(sys.argv[1])
dm2 = skbio.DistanceMatrix.read(sys.argv[2])

if not np.all(dm1.ids == dm2.ids):
    sys.stderr.write("Matrix ids not equal\n")
    sys.exit(1)

if not np.allclose(dm1.data, dm2.data):
    sys.stderr.write("Matrix data not equal\n")
    sys.exit(1)



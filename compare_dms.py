import qiime2
import sys
import skbio


a = qiime2.Artifact.load(sys.argv[1]).view(skbio.DistanceMatrix)
b = qiime2.Artifact.load(sys.argv[2]).view(skbio.DistanceMatrix)

assert a == b

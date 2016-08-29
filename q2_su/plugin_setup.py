# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Emperor development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join

import qiime
import biom
import bp
from q2_types import FeatureTable, Frequency, Phylogeny
from qiime.plugin import Plugin, Properties

import q2_su


def unweighted_unifrac(table : biom.Table,
                       phylogeny : bp.BP) -> skbio.DistanceMatrix:
    pass

plugin = Plugin(
    name='q2_su',
    version=q2_su.__version__,
    website='https://github.com/wasade/q2-state-unifrac',
    package='q2_su',
    citation_text='In preparation',
    user_support_text="Please create an issue on the website tracker"
)

plugin.methods.register_function(
    function=unweighted_unifrac,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny},
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='unweighted UniFrac',
    description="Computes unweighted UniFrac on a user specified table"
)

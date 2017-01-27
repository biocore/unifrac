# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Properties)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.tree import Phylogeny, Rooted

import q2_state_unifrac


plugin = Plugin(
    name='state-unifrac',
    version=q2_state_unifrac.__version__,
    website='https://github.com/wasade/q2-state-unifrac',
    package='q2_state_unifrac',
    user_support_text='https://github.com/wasade/q2-state-unifrac/issues',
    citation_text=None
)


plugin.methods.register_function(
    function=q2_state_unifrac.unweighted,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny[Rooted]},
    parameters={},
    input_descriptions={
        'table': 'A rarefied FeatureTable',
        'phylogeny': ('A rooted phylogeny which relates the observations in '
                      'the table.')
    },
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='Unweighted UniFrac',
    description=('This method computes Unweighted UniFrac')
)


plugin.methods.register_function(
    function=q2_state_unifrac.weighted_unnormalized,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny[Rooted]},
    parameters={},
    input_descriptions={
        'table': 'A rarefied FeatureTable',
        'phylogeny': ('A rooted phylogeny which relates the observations in '
                      'the table.')
    },
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='Weighted unnormalizd UniFrac',
    description=('This method computes weighted unnormalized UniFrac')
)


plugin.methods.register_function(
    function=q2_state_unifrac.weighted_normalized,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny[Rooted]},
    parameters={},
    input_descriptions={
        'table': 'A rarefied FeatureTable',
        'phylogeny': ('A rooted phylogeny which relates the observations in '
                      'the table.')
    },
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='Wweighted normalized UniFrac',
    description=('This method computes weighted normalized UniFrac')
)

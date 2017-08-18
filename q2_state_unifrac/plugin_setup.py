# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Properties, Int, Float, Bool)
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
    parameters={'threads': Int,
                'variance_adjusted': Bool},
    parameter_descriptions={'threads': 'The number of threads to use.',
                            'variance_adjusted':
                                ('Perform variance adjustment based on '
                                 'Chang et al. BMC Bioinformatics 2011')},
    input_descriptions={
        'table': 'A rarefied FeatureTable',
        'phylogeny': ('A rooted phylogeny which relates the observations in '
                      'the table.')
    },
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='Unweighted UniFrac',
    description=('This method computes unweighted UniFrac as described in '
                 'Lozupone and Knight 2005 Appl Environ Microbiol; '
                 'DOI: 10.1128/AEM.71.12.8228-8235.2005')
)


plugin.methods.register_function(
    function=q2_state_unifrac.weighted_unnormalized,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny[Rooted]},
    parameters={'threads': Int,
                'variance_adjusted': Bool},
    parameter_descriptions={'threads': 'The number of threads to use.',
                            'variance_adjusted':
                                ('Perform variance adjustment based on '
                                 'Chang et al. BMC Bioinformatics 2011')},
    input_descriptions={
        'table': 'A rarefied FeatureTable',
        'phylogeny': ('A rooted phylogeny which relates the observations in '
                      'the table.')
    },
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='Weighted unnormalizd UniFrac',
    description=('This method computes weighted unnormalized UniFrac as '
                 'described in Lozupone et al. 2007 Appl Environ Microbiol; '
                 'DOI: 10.1128/AEM.01996-06')
)


plugin.methods.register_function(
    function=q2_state_unifrac.weighted_normalized,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny[Rooted]},
    parameters={'threads': Int,
                'variance_adjusted': Bool},
    parameter_descriptions={'threads': 'The number of threads to use.',
                            'variance_adjusted':
                                ('Perform variance adjustment based on '
                                 'Chang et al. BMC Bioinformatics 2011')},
    input_descriptions={
        'table': 'A rarefied FeatureTable',
        'phylogeny': ('A rooted phylogeny which relates the observations in '
                      'the table.')
    },
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='Weighted normalized UniFrac',
    description=('This method computes weighted normalized UniFrac as '
                 'described in Lozupone et al. 2007 Appl Environ Microbiol; '
                 'DOI: 10.1128/AEM.01996-06')
)


plugin.methods.register_function(
    function=q2_state_unifrac.generalized,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny[Rooted]},
    parameters={'threads': Int,
                'variance_adjusted': Bool,
                'alpha': Float},
    parameter_descriptions={'threads': 'The number of threads to use.',
                            'variance_adjusted':
                                ('Perform variance adjustment based on '
                                 'Chang et al. BMC Bioinformatics 2011'),
                            'alpha': ('The value of alpha controls importance '
                                      'of sample proportions. 1.0 is '
                                      'weighted normalized UniFrac. 0.0 is '
                                      'close to unweighted UniFrac, but only '
                                      'if the sample proportions are '
                                      'dichotomized.')},
    input_descriptions={
        'table': 'A rarefied FeatureTable',
        'phylogeny': ('A rooted phylogeny which relates the observations in '
                      'the table.')
    },
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='Generalized UniFrac',
    description=('This method computes Generalized UniFrac as described in '
                 'Chen et al. 2012 Bioinformatics; '
                 'DOI: 10.1093/bioinformatics/bts342')
)

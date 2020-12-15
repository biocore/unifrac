# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, UniFrac development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from unifrac._methods import (unweighted,
                              weighted_normalized,
                              weighted_unnormalized,
                              generalized,
                              unweighted_fp32,
                              weighted_normalized_fp32,
                              weighted_unnormalized_fp32,
                              generalized_fp32, 
                              meta)
from unifrac._api import ssu, faith_pd, ssu_to_file


__version__ = pkg_resources.get_distribution('unifrac').version
__all__ = ['unweighted', 'weighted_normalized', 'weighted_unnormalized',
           'generalized', 'unweighted_p32', 'weighted_normalized_fp32',
           'weighted_unnormalized_fp32', 'generalized_fp32',
           'meta', 
           'ssu', 'faith_pd', 'ssu_to_file']

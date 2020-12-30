# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, UniFrac development team.
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
                              unweighted_to_file,
                              weighted_normalized_to_file,
                              weighted_unnormalized_to_file,
                              generalized_to_file,
                              unweighted_fp32_to_file,
                              weighted_normalized_fp32_to_file,
                              weighted_unnormalized_fp32_to_file,
                              generalized_fp32_to_file,
                              meta,
                              h5unifrac,
                              h5pcoa)
from unifrac._api import ssu, faith_pd, ssu_to_file


__version__ = pkg_resources.get_distribution('unifrac').version
__all__ = ['unweighted', 'weighted_normalized', 'weighted_unnormalized',
           'generalized', 'unweighted_fp32', 'weighted_normalized_fp32',
           'weighted_unnormalized_fp32', 'generalized_fp32',
           'meta',
           'unweighted_to_file', 'weighted_normalized_to_file',
           'weighted_unnormalized_to_file',
           'generalized_to_file', 'unweighted_fp32_to_file',
           'weighted_normalized_fp32_to_file',
           'weighted_unnormalized_fp32_to_file',
           'generalized_fp32_to_file',
           'h5unifrac', 'h5pcoa',
           'ssu', 'faith_pd', 'ssu_to_file']

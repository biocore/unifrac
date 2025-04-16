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
                              unweighted_fp64,
                              weighted_normalized_fp64,
                              weighted_unnormalized_fp64,
                              generalized_fp64,
                              unweighted_fp32,
                              weighted_normalized_fp32,
                              weighted_unnormalized_fp32,
                              generalized_fp32,
                              unweighted_to_file,
                              weighted_normalized_to_file,
                              weighted_unnormalized_to_file,
                              generalized_to_file,
                              unweighted_fp64_to_file,
                              weighted_normalized_fp64_to_file,
                              weighted_unnormalized_fp64_to_file,
                              generalized_fp64_to_file,
                              unweighted_fp32_to_file,
                              weighted_normalized_fp32_to_file,
                              weighted_unnormalized_fp32_to_file,
                              generalized_fp32_to_file,
                              meta,
                              h5unifrac, h5unifrac_all,
                              h5pcoa, h5pcoa_all,
                              h5permanova, h5permanova_dict,
                              faith_pd)
from unifrac._api import ssu, ssu_fast, set_random_seed
from unifrac._api import ssu_to_file, ssu_to_file_v2, ssu_inmem
from unifrac._api import faith_pd as _faith_pd  # noqa: F401

__version__ = pkg_resources.get_distribution('unifrac').version
__all__ = ['unweighted', 'weighted_normalized', 'weighted_unnormalized',
           'generalized', 'unweighted_fp64', 'weighted_normalized_fp64',
           'weighted_unnormalized_fp64', 'generalized_fp64',
           'unweighted_fp32', 'weighted_normalized_fp32',
           'weighted_unnormalized_fp32', 'generalized_fp32',
           'meta', 'set_random_seed',
           'unweighted_to_file', 'weighted_normalized_to_file',
           'weighted_unnormalized_to_file',
           'generalized_to_file', 'unweighted_fp64_to_file',
           'weighted_normalized_fp64_to_file',
           'weighted_unnormalized_fp64_to_file',
           'generalized_fp64_to_file',
           'unweighted_fp32_to_file',
           'weighted_normalized_fp32_to_file',
           'weighted_unnormalized_fp32_to_file',
           'generalized_fp32_to_file',
           'h5unifrac', 'h5unifrac_all', 'h5pcoa', 'h5pcoa_all',
           'h5permanova', 'h5permanova_dict',
           'ssu', 'ssu_fast', 'faith_pd',
           'ssu_to_file', 'ssu_to_file_v2',
           'ssu_inmem']

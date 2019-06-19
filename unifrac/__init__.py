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
                              generalized, meta)
from unifrac._api import ssu, faith_pd


__version__ = pkg_resources.get_distribution('unifrac').version
__all__ = ['unweighted', 'weighted_normalized', 'weighted_unnormalized',
           'generalized', 'meta', 'ssu', 'faith_pd']

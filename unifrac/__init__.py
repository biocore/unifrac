# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from q2_state_unifrac._methods import (unweighted,
                                       weighted_normalized,
                                       weighted_unnormalized,
                                       generalized, meta)
from q2_state_unifrac._api import ssu


__version__ = pkg_resources.get_distribution('q2-state-unifrac').version
__all__ = ['unweighted', 'weighted_normalized', 'weighted_unnormalized',
           'generalized', 'meta', 'ssu']

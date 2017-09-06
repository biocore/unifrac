# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, UniFrac development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# Code pulled from cogent.maths.unifrac.fast_unifrac; the authors have
# previously indicated approval for converstion from GPL -> BSD
# https://github.com/biocore/scikit-bio#the-pre-history-of-scikit-bio

# These methods did not have unit tests in cogent

import numpy as np


def consolidate_skipping_missing_matrices(matrices, env_names, weights,
                                          all_env_names):
    """Consolidates matrices, skipping any that are missing envs"""
    weight_sum = 0
    result = np.zeros((len(all_env_names), len(all_env_names)), float)
    for m, e, w in zip(matrices, env_names, weights):
        if e == all_env_names:  # note -- assumes sorted
            result += m * w
            weight_sum += w
    # readjust weights for missing matrices
    result /= weight_sum
    return result


def consolidate_missing_zero(matrices, env_names, weights, all_env_names):
    """Consolidates matrices, setting missing values to 0 distance"""
    result = np.zeros((len(all_env_names), len(all_env_names)), float)
    for m, e, w in zip(matrices, env_names, weights):
        result += reshape_by_name(m, e, all_env_names, 0) * w
    return result


def consolidate_missing_one(matrices, env_names, weights, all_env_names):
    """Consolidates matrices, setting missing values to 1 distance"""
    result = np.zeros((len(all_env_names), len(all_env_names)), float)
    for m, e, w in zip(matrices, env_names, weights):
        result += reshape_by_name(m, e, all_env_names, 1) * w
    return result


def consolidate_skipping_missing_values(matrices, env_names, weights,
                                        all_env_names):
    """Consolidates matrices, skipping only values from missing envs"""
    result = []
    for m, e, w in zip(matrices, env_names, weights):
        reshaped = reshape_by_name(m, e, all_env_names, masked=True)
        reshaped *= w
        result.append(reshaped)

    data = np.array([i.data for i in result], float)
    masks = np.array([i.mask for i in result], bool)
    masked_result = np.ma.array(data, mask=masks)
    # figure out mask of weights so we can figure out per-element weighting
    masked_weights = np.ma.array(np.zeros(data.shape), mask=masks) + \
        np.array(weights, float).reshape((len(weights), 1, 1))
    return masked_result.sum(0) / masked_weights.sum(0)


def reshape_by_name(m, old_names, new_names, default_off_diag=0,
                    default_diag=0, masked=False):
    """Reshape matrix m mapping slots from old names to new names. """
    num_names = len(new_names)
    result = np.zeros((num_names, num_names), float) + default_off_diag
    for i in range(num_names):
        result[i, i] = default_diag
    pairs = {}
    for i, n in enumerate(old_names):
        if n in new_names:
            pairs[i] = new_names.index(n)
    for i, row in enumerate(m):
        new_i = pairs[i]
        for j, val in enumerate(row):
            new_j = pairs[j]
            result[new_i, new_j] = val
    if masked:
        mask = np.ones((num_names, num_names), float)
        for i in pairs.values():
            for j in pairs.values():
                mask[i, j] = 0
        result = np.ma.array(result, mask=mask)
    return result


CONSOLIDATIONS = \
    {'skipping_missing_matrices': consolidate_skipping_missing_matrices,
     'missing_zero': consolidate_missing_zero,
     'missing_one': consolidate_missing_one,
     'skipping_missing_values': consolidate_skipping_missing_values}

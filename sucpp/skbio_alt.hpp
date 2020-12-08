/*
 * Classes, methods and unction that provide skbio-like unctionality
 */

#ifndef UNIFRAC_SKBIO_ALT_H
#define UNIFRAC_SKBIO_ALT_H

#include <stdint.h>

namespace su {

// Center the matrix
// mat and center must be nxn and symmetric
// centered must be pre-allocated and same size as mat...will work even if centered==mat
void mat_to_centered(const double * mat, const uint32_t n_samples, double * centered);
void mat_to_centered(const float  * mat, const uint32_t n_samples, float  * centered);
void mat_to_centered(const double * mat, const uint32_t n_samples, float  * centered);

// Find eigen values and vectors
// Based on N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
//     Original Paper: https://arxiv.org/abs/1007.5510
// centered == n x n, must be symmetric, Note: will be used in-place as temp buffer
void find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, double * centered, double * &eigenvalues, double * &eigenvectors);
void find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, float  * centered, float  * &eigenvalues, float  * &eigenvectors);

// Perform Principal Coordinate Analysis
// mat       - in, result of unifrac compute
// n_samples - in, size of the matrix (n x n)
// n_dims    - in, Dimensions to reduce the distance matrix to. This number determines how many eigenvectors and eigenvalues will be returned.
// eigenvalues - out, alocated buffer of size n_dims
// samples     - out, alocated buffer of size n_dims x n_samples
// proportion_explained - out, allocated buffer of size n_dims
void pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, double * &eigenvalues, double * &samples, double * &proportion_explained);
void pcoa(const float  * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained);
void pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained);

// in-place version, will use mat as temp buffer internally
void pcoa_inplace(double * mat, const uint32_t n_samples, const uint32_t n_dims, double * &eigenvalues, double * &samples, double * &proportion_explained);
void pcoa_inplace(float  * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained);


}

#endif

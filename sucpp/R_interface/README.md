# R interface for Strided State Unifrac

This provides an R interface for Unweighted Unifrac. This interface works using
R's Rcpp library. To load this, in R use `library(Rcpp)` and
`sourceCpp("su_R.cpp")`. The Unifrac method takes in three arguments: a file
path to an HDF5 formatted BIOM table, a filepath to a newick formatted tree
file, and the number of threads to be used It is expected that the observations
described in the BIOM table correspond to a subset of the tips of the input
tree. The method returns a list containing an `int` `n_samples`, denoting the
number of samples in the table, a `boolean` `is_upper_triangle`, denoting 
whether Unifrac generated a square matrix and if it has returned the upper 
triangle , an `int` `cf_size`, denoting the size of the condensed form of the 
matrix, and `c_form`, an array representation of the condensed form of the 
matrix, obtained by taking the upper triangle.

```R
> library(Rcpp)
> sourceCpp("su_R.cpp")
> table = "../test.biom"
> tree = "../test.tre"
> nthreads = 2
> unif = unifrac(table, tree, nthreads)
> unif
$n_samples
[1] 6

$is_sqaure
[1] TRUE

$cf_size
[1] 15

$c_form
 [1] 0.2000000 0.5714286 0.6000000 0.5000000 0.2000000 0.4285714 0.6666667
 [8] 0.6000000 0.3333333 0.7142857 0.8571429 0.4285714 0.3333333 0.4000000
[15] 0.6000000

```
 

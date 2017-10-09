# R interface for Strided State Unifrac

This provides an R interface for Unweighted Unifrac. This interface works using
R's Rcpp library. To load this, in R use `library(Rcpp)` and 
`sourceCpp("su_R.cpp")`. The unifrac method takes in three arguments, a biom
table, and tree and the number of threads to be used. The method returns a list
containing an `int` `n_samples`, denoting the number of samples in the table, 
a `boolean` `is_sqaure`, denoting whether Unifrac generated a square matrix, 
an `int` `cf_size`, denoting the size of the condensed form of the matrix, 
and `c_form`, an array representation of the condensed form of the matrix,
obtained by taking the upper triangle. 

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
 

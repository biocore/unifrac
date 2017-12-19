library(Rcpp)

equals <- function(x, y){
	stopifnot(x==y)
}
aboutEquals <- function(x, y){
	stopifnot((x-y)<0.005)
}
source = "su_R.cpp"
sourceCpp(source)
table = "test.biom"
tree = "test.tre"
nthreads = 2
unif = unifrac(table, tree, nthreads)

exp = c(0.2000000, 0.5714286, 0.6000000, 0.5000000, 0.2000000, 
		0.4285714, 0.6666667, 0.6000000, 0.3333333, 0.7142857, 
		0.8571429, 0.4285714, 0.3333333, 0.4000000, 0.6000000)


equals(unif["n_samples"][[1]], 6)
equals(unif["cf_size"][[1]], 15)
equals(unif["is_upper_triangle"][[1]], TRUE)


for ( i in 1:15){
	aboutEquals(unif["c_form"][[1]][i], exp[i])
}
print('All tests pass')

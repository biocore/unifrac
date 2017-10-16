library(Rcpp)
library(RUnit)
sourceCpp("su_R.cpp")
table = "../test.biom"
tree = "../test.tre"
nthreads = 2
unif = unifrac(table, tree, nthreads)

exp = c(0.2000000, 0.5714286, 0.6000000, 0.5000000, 0.2000000, 
		0.4285714, 0.6666667, 0.6000000, 0.3333333, 0.7142857, 
		0.8571429, 0.4285714, 0.3333333, 0.4000000, 0.6000000)


checkEquals(unif["n_samples"][[1]], 6)
checkEquals(unif["cf_size"][[1]], 15)
checkTrue(unif["is_sqaure"][[1]])

print("testing values of c_form")
for ( i in 1:15){
	print(checkEquals(unif["c_form"][[1]][i], exp[i], tolerance=0.005))
}
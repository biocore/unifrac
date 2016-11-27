set -e
set -x

basedir=bench_tables_trees
resdir=$basedir/results
mkdir -p $resdir
for f in $basedir/*.biom
do
    bench=${basedir}/$(basename $f .biom)
    res=${resdir}/$(basename $f .biom)
    for method in {unweighted,weighted_normalized,weighted_unnormalized}
    do
        /usr/bin/time -l ./su ${bench}.tre ${bench}.biom $method > ${res}.${method}.su.dm 2> ${res}.${method}.su.stats
        /usr/bin/time -l ./sk ${bench}.tre ${bench}.biom $method > ${res}.${method}.sk.dm 2> ${res}.${method}.sk.stats
        python compare_dms.py ${res}.${method}.sk.dm ${res}.${method}.su.dm
    done
done

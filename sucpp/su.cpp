#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include "cmd.hpp"
#include <thread>

void usage() {
    std::cout << "usage: ssu -i <biom> -o <out.dm> -m [METHOD] -t <newick> [-n threads]" << std::endl;
    std::cout << std::endl;
    std::cout << "    -i\tThe input BIOM table." << std::endl;
    std::cout << "    -t\tThe input phylogeny in newick." << std::endl;
    std::cout << "    -m\tThe method, [unweighted | weighted_normalized | weighted_unnormalized]." << std::endl;
    std::cout << "    -o\tThe output distance matrix." << std::endl;
    std::cout << "    -n\t[OPTIONAL] The number of threads, default is 1." << std::endl;
    std::cout << std::endl;
}

void err(std::string msg) {
    std::cerr << "ERROR: " << msg << std::endl << std::endl;
    usage();
}

int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        usage();
        return EXIT_SUCCESS;
    }

    if(!input.cmdOptionExists("-i")) {
        usage();
        return EXIT_SUCCESS;
    } 
    if(!input.cmdOptionExists("-t")) {
        usage();
        return EXIT_SUCCESS;
    } 
    if(!input.cmdOptionExists("-o")) {
        usage();
        return EXIT_SUCCESS;
    } 
    if(!input.cmdOptionExists("-m")) {
        usage();
        return EXIT_SUCCESS;
    } 

    unsigned int nthreads;
    const std::string &table_filename = input.getCmdOption("-i");
    const std::string &tree_filename = input.getCmdOption("-t");
    const std::string &output_filename = input.getCmdOption("-o");
    const std::string &method_string = input.getCmdOption("-m");
    const std::string &nthreads_arg = input.getCmdOption("-n");

    su::Method method;

    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }

    if(tree_filename.empty()) {
        err("tree filename missing");
        return EXIT_FAILURE;
    }

    if(output_filename.empty()) {
        err("output filename missing");
        return EXIT_FAILURE;
    }
    if(method_string.empty()) {
        err("method missing");
        return EXIT_FAILURE;
    }
    if(nthreads_arg.empty()) {
        nthreads = 1;
    } else {
        nthreads = atoi(nthreads_arg.c_str());
    }
    
    if(method_string == "unweighted") 
        method = su::unweighted;
    else if(method_string == "weighted_normalized")
        method = su::weighted_normalized;
    else if(method_string == "weighted_unnormalized")
        method = su::weighted_unnormalized;
    else {
        err("Unknown method");
        return EXIT_FAILURE;
    }

    std::ifstream ifs((char*)(tree_filename).c_str());
    std::string content = std::string(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
    su::biom table = su::biom(table_filename);
    su::BPTree tree = su::BPTree(content);
    
    std::unordered_set<std::string> to_keep(table.obs_ids.begin(), table.obs_ids.end());
    su::BPTree tree_sheared = tree.shear(to_keep).collapse();

    std::vector<double*> dm_stripes; // = su::make_strides(table.n_samples);
    std::vector<double*> dm_stripes_total; // = su::make_strides(table.n_samples);
    dm_stripes.resize((table.n_samples + 1) / 2);
    dm_stripes_total.resize((table.n_samples + 1) / 2);
    
    unsigned int chunksize = dm_stripes.size() / nthreads;
    unsigned int start = 0;
    unsigned int end = dm_stripes.size();
    unsigned int *starts = (unsigned int*)malloc(sizeof(unsigned int) * nthreads);
    unsigned int *ends = (unsigned int*)malloc(sizeof(unsigned int) * nthreads);
    std::vector<std::thread> threads(nthreads);
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        starts[tid] = start;
        ends[tid] = start + chunksize;
        start = start + chunksize;
    }
    // the last thread gets any trailing bits
    ends[threads.size() - 1] = end;
    
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        threads[tid] = std::thread(su::unifrac, 
                                   std::ref(table),
                                   std::ref(tree_sheared), 
                                   method, 
                                   std::ref(dm_stripes), 
                                   std::ref(dm_stripes_total), 
                                   starts[tid], 
                                   ends[tid], 
                                   tid);
    }

    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        threads[tid].join();
    }
    
    free(starts);
    free(ends);

    double **dm = su::deconvolute_stripes(dm_stripes, table.n_samples);

    std::ofstream output;
    output.open(output_filename);
    
    for(unsigned int i = 0; i < table.n_samples; i++)
        output << "\t" << table.sample_ids[i];
    output << std::endl;
    for(unsigned int i = 0; i < table.n_samples; i++) {
        output << table.sample_ids[i];
        for(unsigned int j = 0; j < table.n_samples; j++)
            output << std::setprecision(16) << "\t" << dm[i][j];
        output << std::endl;
    }
    unsigned int n_rotations = (table.n_samples + 1) / 2;
    for(unsigned int i = 0; i < n_rotations; i++)
        free(dm_stripes[i]);
    
    if(method == su::weighted_normalized || method == su::unweighted) {
        for(unsigned int i = 0; i < n_rotations; i++)
            free(dm_stripes_total[i]);
    }
    return EXIT_SUCCESS;
}

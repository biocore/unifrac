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
    std::cout << "ya, you know, make this" << std::endl;
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

    std::vector<std::thread> threads(nthreads);
    double **dm = unifrac(table, tree_sheared, method, threads);

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
    return EXIT_SUCCESS;
}

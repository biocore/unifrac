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
    std::cout << "usage: ssu -i <biom> -o <out.dm> -m [METHOD] -t <newick> [-n threads] [-a alpha] [--vaw]" << std::endl;
    std::cout << std::endl;
    std::cout << "    -i\t\tThe input BIOM table." << std::endl;
    std::cout << "    -t\t\tThe input phylogeny in newick." << std::endl;
    std::cout << "    -m\t\tThe method, [unweighted | weighted_normalized | weighted_unnormalized | generalized]." << std::endl;
    std::cout << "    -o\t\tThe output distance matrix." << std::endl;
    std::cout << "    -n\t\t[OPTIONAL] The number of threads, default is 1." << std::endl;
    std::cout << "    -a\t\t[OPTIONAL] Generalized UniFrac alpha, default is 1." << std::endl;
    std::cout << "    --vaw\t[OPTIONAL] Variance adjusted, default is to not adjust for variance." << std::endl;
    std::cout << std::endl;
    std::cout << "Citations: " << std::endl;
    std::cout << "    For UniFrac, please see:" << std::endl;
    std::cout << "        Lozupone and Knight Appl Environ Microbiol 2005; DOI: 10.1128/AEM.71.12.8228-8235.2005" << std::endl;
    std::cout << "        Lozupone et al. Appl Environ Microbiol 2007; DOI: 10.1128/AEM.01996-06" << std::endl;
    std::cout << "        Hamady et al. ISME 2010; DOI: 10.1038/ismej.2009.97" << std::endl;
    std::cout << "        Lozupone et al. ISME 2011; DOI: 10.1038/ismej.2010.133" << std::endl;
    std::cout << "    For Generalized UniFrac, please see: " << std::endl;
    std::cout << "        Chen et al. Bioinformatics 2012; DOI: 10.1093/bioinformatics/bts342" << std::endl;
    std::cout << "    For Variance Adjusted UniFrac, please see: " << std::endl;
    std::cout << "        Chang et al. BMC Bioinformatics 2011; DOI: 10.1186/1471-2105-12-118" << std::endl;
    std::cout << std::endl;
}

void err(std::string msg) {
    std::cerr << "ERROR: " << msg << std::endl << std::endl;
    usage();
}

// https://stackoverflow.com/a/19841704/19741
bool is_file_exists(const char *fileName) {
    std::ifstream infile(fileName);
        return infile.good();
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
    const std::string &gunifrac_arg = input.getCmdOption("-a");

    su::Method method;

    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }
    if(!is_file_exists(table_filename.c_str())) {
        return EXIT_FAILURE;
    }

    if(tree_filename.empty()) {
        err("tree filename missing");
        return EXIT_FAILURE;
    }
    if(!is_file_exists(tree_filename.c_str())) {
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
    
    bool vaw = input.cmdOptionExists("--vaw"); 
    double g_unifrac_alpha;
    if(gunifrac_arg.empty()) {
        g_unifrac_alpha = 1.0;
    } else {
        g_unifrac_alpha = atof(gunifrac_arg.c_str());
    }

    if(method_string == "unweighted") 
        method = su::unweighted;
    else if(method_string == "weighted_normalized")
        method = su::weighted_normalized;
    else if(method_string == "weighted_unnormalized")
        method = su::weighted_unnormalized;
    else if(method_string == "generalized")
        method = su::generalized;
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

    std::vector<double*> dm_stripes; 
    std::vector<double*> dm_stripes_total; 
    dm_stripes.resize((table.n_samples + 1) / 2);
    dm_stripes_total.resize((table.n_samples + 1) / 2);
    
    unsigned int chunksize = dm_stripes.size() / nthreads;
    unsigned int start = 0;
    unsigned int end = dm_stripes.size();
    
    std::vector<su::task_parameters> tasks;
    tasks.resize(nthreads);

    std::vector<std::thread> threads(nthreads);
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        tasks[tid].tid = tid;
        tasks[tid].start = start; // stripe start
        tasks[tid].stop = start + chunksize;  // stripe end 
        tasks[tid].n_samples = table.n_samples;
        tasks[tid].g_unifrac_alpha = g_unifrac_alpha;
        start = start + chunksize;
    }
    // the last thread gets any trailing bits
    tasks[threads.size() - 1].stop = end;
    
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        if(vaw)
            threads[tid] = std::thread(su::unifrac_vaw, 
                                       std::ref(table),
                                       std::ref(tree_sheared), 
                                       method, 
                                       std::ref(dm_stripes), 
                                       std::ref(dm_stripes_total), 
                                       &tasks[tid]);
        else
            threads[tid] = std::thread(su::unifrac, 
                                       std::ref(table),
                                       std::ref(tree_sheared), 
                                       method, 
                                       std::ref(dm_stripes), 
                                       std::ref(dm_stripes_total), 
                                       &tasks[tid]);
    }

    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        threads[tid].join();
    }
    
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

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "api.hpp"
#include "cmd.hpp"
#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"

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
    
    bool vaw = input.cmdOptionExists("--vaw"); 
    double g_unifrac_alpha;
    if(gunifrac_arg.empty()) {
        g_unifrac_alpha = 1.0;
    } else {
        g_unifrac_alpha = atof(gunifrac_arg.c_str());
    }

    mat_t *result = NULL;
    compute_status status;
    status = one_off(table_filename.c_str(), tree_filename.c_str(), method_string.c_str(), 
                     vaw, g_unifrac_alpha, nthreads, &result);
    if(status != okay || result == NULL) {
        fprintf(stderr, "Compute failed in one_off with error code: %d\n", status);
        exit(EXIT_FAILURE);
    }
    
    std::ofstream output;
    output.open(output_filename);
  
    uint32_t comb_N = su::comb_2(result->n_samples);
    uint32_t comb_N_minus;
    double v;
    for(unsigned int i = 0; i < result->n_samples; i++)
        output << "\t" << result->sample_ids[i];
    output << std::endl;
    for(unsigned int i = 0; i < result->n_samples; i++) {
        output << result->sample_ids[i];
        for(unsigned int j = 0; j < result->n_samples; j++) {
            if(i < j) { // upper triangle
                comb_N_minus = su::comb_2(result->n_samples - i);
                v = result->condensed_form[comb_N - comb_N_minus + (j - i - 1)];
            } else if (i > j) { // lower triangle
                comb_N_minus = su::comb_2(result->n_samples - j);
                v = result->condensed_form[comb_N - comb_N_minus + (i - j - 1)];
            } else {
                v = 0.0;
            }
            output << std::setprecision(16) << "\t" << v;
        }
        output << std::endl;
    }
    destroy_mat(&result);

    return EXIT_SUCCESS;
}

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
    std::cout << "    [--mode [MODE]] [--start starting-stripe] [--stop stopping-stripe]" << std::endl;
    std::cout << std::endl;
    std::cout << "    -i\t\tThe input BIOM table." << std::endl;
    std::cout << "    -t\t\tThe input phylogeny in newick." << std::endl;
    std::cout << "    -m\t\tThe method, [unweighted | weighted_normalized | weighted_unnormalized | generalized]." << std::endl;
    std::cout << "    -o\t\tThe output distance matrix." << std::endl;
    std::cout << "    -n\t\t[OPTIONAL] The number of threads, default is 1." << std::endl;
    std::cout << "    -a\t\t[OPTIONAL] Generalized UniFrac alpha, default is 1." << std::endl;
    std::cout << "    --vaw\t[OPTIONAL] Variance adjusted, default is to not adjust for variance." << std::endl;
    std::cout << "    --mode\t[OPTIONAL] Mode of operation:" << std::endl;
    std::cout << "    \t\t    one-off : [DEFAULT] compute UniFrac." << std::endl;
    std::cout << "    \t\t    partial : Compute UniFrac over a subset of stripes." << std::endl;
    std::cout << "    \t\t    merge-partial : Merge partial UniFrac results." << std::endl;
    std::cout << "    \t\t    stripes : Print the number of stripes." << std::endl;
    std::cout << "    --start\t[OPTIONAL] If mode==partial, the starting stripe." << std::endl;
    std::cout << "    --stop\t[OPTIONAL] If mode==partial, the stopping stripe." << std::endl;
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

int mode_partial(std::string table_filename, std::string tree_filename, 
                 std::string output_filename, std::string method_string,
                 bool vaw, double g_unifrac_alpha, unsigned int nthreads,
                 int start_stripe, int stop_stripe) {
    if(start_stripe < 0)
        err("Starting stripe must be >= 0");
    if(stop_stripe <= start_stripe)
        err("Stopping stripe must be > start stripe");

    partial_mat_t *result = NULL;
    compute_status status;
    status = partial(table_filename.c_str(), tree_filename.c_str(), method_string.c_str(), 
                     vaw, g_unifrac_alpha, nthreads, start_stripe, stop_stripe, &result);
    if(status != okay || result == NULL) {
        fprintf(stderr, "Compute failed in one_off with error code: %d\n", status);
        exit(EXIT_FAILURE);
    }
    
    std::ofstream output;
    output.open(output_filename, std::ios::binary);
 
    std::string magic("SSU-PARTIAL-00");
    output << magic;
    output.write(reinterpret_cast<const char*>(&result->n_samples), sizeof(uint32_t));
    // // // // /// // WRITE OUT what stripes are present up front, so a collection of files 
    // can be surveyed rapidly to test for consistency
    for(unsigned int i = 0; i < result->n_samples; i++) {
        uint16_t length = strlen(result->sample_ids[i]);
        output.write(reinterpret_cast<const char*>(&length), sizeof(uint16_t));
        output << result->sample_ids[i];
    }
    uint32_t n_stripes = result->stripe_stop - result->stripe_start;
    output.write(reinterpret_cast<const char*>(&n_stripes), sizeof(uint32_t));
    for(unsigned int i = 0; i < n_stripes; i++) {
        uint32_t current_stripe = i + result->stripe_start;
        output.write(reinterpret_cast<const char*>(&current_stripe), sizeof(uint32_t));
        
        /// :( streamsize didn't seem to work. probably a fancy way to do this, but the regular loop is fast too
        //output.write(reinterpret_cast<const char*>(&result->stripes[i]), std::streamsize(sizeof(double) * result->n_samples));
        for(unsigned int j = 0; j < result->n_samples; j++)
            output.write(reinterpret_cast<const char*>(&result->stripes[i][j]), sizeof(double));
    }
    // // // // // WRITE OUT an end
    destroy_partial_mat(&result);

    return EXIT_SUCCESS;
}

int mode_one_off(std::string table_filename, std::string tree_filename, 
                 std::string output_filename, std::string method_string,
                 bool vaw, double g_unifrac_alpha, unsigned int nthreads) {

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
    const std::string &mode_arg = input.getCmdOption("--mode");
    const std::string &start_arg = input.getCmdOption("--start");
    const std::string &stop_arg = input.getCmdOption("--stop");

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

    int start_stripe;
    if(start_arg.empty()) 
        start_stripe = 0;
    else
        start_stripe = atoi(start_arg.c_str());

    int stop_stripe;
    if(stop_arg.empty()) 
        stop_stripe = 0;
    else
        stop_stripe = atoi(stop_arg.c_str());

    if(mode_arg.empty() || mode_arg == "one-off")
        return mode_one_off(table_filename, tree_filename, output_filename, method_string, vaw, g_unifrac_alpha, nthreads);
    else if(mode_arg == "partial")
        return mode_partial(table_filename, tree_filename, output_filename, method_string, vaw, g_unifrac_alpha, nthreads, start_stripe, stop_stripe);
    else if(mode_arg == "merge-partial")
        err("Not implemented");
    else if(mode_arg == "stripes")
        err("Not implemented");
    else 
        err("Unknown mode. Valid options are: one-off, partial, merge-partial");

    return EXIT_SUCCESS;
}


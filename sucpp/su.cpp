#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <glob.h>
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
    std::cout << "    -f\t\t[OPTIONAL] Bypass tips, reduces compute by about 50%." << std::endl;
    std::cout << "    --vaw\t[OPTIONAL] Variance adjusted, default is to not adjust for variance." << std::endl;
    std::cout << "    --mode\t[OPTIONAL] Mode of operation:" << std::endl;
    std::cout << "    \t\t    one-off : [DEFAULT] compute UniFrac." << std::endl;
    std::cout << "    \t\t    partial : Compute UniFrac over a subset of stripes." << std::endl;
    std::cout << "    \t\t    partial-report : Start and stop suggestions for partial compute." << std::endl;
    std::cout << "    \t\t    merge-partial : Merge partial UniFrac results." << std::endl;
    std::cout << "    --start\t[OPTIONAL] If mode==partial, the starting stripe." << std::endl;
    std::cout << "    --stop\t[OPTIONAL] If mode==partial, the stopping stripe." << std::endl;
    std::cout << "    --partial-pattern\t[OPTIONAL] If mode==merge-partial, a glob pattern for partial outputs to merge." << std::endl;
    std::cout << "    --n-partials\t[OPTIONAL] If mode==partial-report, the number of partitions to compute." << std::endl;
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
    std::cout << "Runtime progress can be obtained by issuing a SIGUSR1 signal. If running with " << std::endl;
    std::cout << "multiple threads, this signal will only be honored if issued to the master PID. " << std::endl;
    std::cout << "The report will yield the following information: " << std::endl;
    std::cout << std::endl;
    std::cout << "tid:<thread ID> start:<starting stripe> stop:<stopping stripe> k:<postorder node index> total:<number of nodes>" << std::endl;
    std::cout << std::endl;
    std::cout << "The proportion of the tree that has been evaluated can be determined from (k / total)." << std::endl;
    std::cout << std::endl;
}


// https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
inline std::vector<std::string> glob(const std::string& pat){
    using namespace std;
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        std::cout << glob_result.gl_pathv[i] << std::endl;
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}


void err(std::string msg) {
    std::cerr << "ERROR: " << msg << std::endl << std::endl;
    usage();
}

int mode_partial_report(const std::string table_filename, int npartials) {
    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }

    su::biom table = su::biom(table_filename.c_str());
    int total_stripes = (table.n_samples + 1) / 2;
    std::cout << "Total samples: " << table.n_samples << std::endl;
    std::cout << "Total stripes: " << total_stripes << std::endl;

    unsigned int fullchunk = (total_stripes + npartials - 1) / npartials;  // this computes the ceiling
    unsigned int smallchunk = total_stripes / npartials;
    
    unsigned int n_fullbins = total_stripes % npartials;
    if(n_fullbins == 0)
        n_fullbins = npartials;

    unsigned int start = 0;
    unsigned int stop = 0;
    for(unsigned int p = 0; p < npartials; p++) {
        if(p < n_fullbins) {
            stop = start + fullchunk;
            std::cout << "Partition " << p << ", suggested start and stop: " << start << ", " << stop << std::endl;
            start = start + fullchunk;
        } else {
            stop = start + smallchunk;  // stripe end 
            std::cout << "Partition " << p << ", suggested start and stop: " << start << ", " << stop << std::endl;
            start = start + smallchunk;
        }
    }
} 

int mode_merge_partial(std::string output_filename,
					   std::string partial_pattern,
					   unsigned int nthreads) {
    if(output_filename.empty()) {
        err("output filename missing");
        return EXIT_FAILURE;
    }

    if(partial_pattern.empty()) {
        std::string msg("Partial file pattern missing. For instance, if your partial results\n" \
                        "are named 'ssu.unweighted.start0.partial', 'ssu.unweighted.start10.partial', \n" \
                        "etc, then a pattern of 'ssu.unweighted.start*.partial' would make sense");
        err(msg);
        return EXIT_FAILURE;
    }
    
    mat_t *result = NULL;
    
    /// glob is not working with lldb. not sure why yet.
    std::vector<std::string> partials = glob(partial_pattern);
    partial_mat_t** partial_mats = (partial_mat_t**)malloc(sizeof(partial_mat_t*) * partials.size());
    for(size_t i = 0; i < partials.size(); i++) {
        IOStatus io_err = read_partial(partials[i].c_str(), &partial_mats[i]);
        if(io_err != read_okay) {
            std::ostringstream msg;
            msg << "Unable to parse file (" << partials[i] << "); err " << io_err;
            err(msg.str());
            return EXIT_FAILURE;
        }
    }

    MergeStatus status = merge_partial(partial_mats, partials.size(), nthreads, &result);
    
    if(status != merge_okay) {
        std::ostringstream msg;
        msg << "Unable to complete merge; err " << status;
        err(msg.str());
        return EXIT_FAILURE;
    }

    IOStatus io_err = write_mat(output_filename.c_str(), result);
    if(io_err != write_okay) {
        std::ostringstream msg;
        msg << "Unable to write; err " << io_err;
        err(msg.str());
        return EXIT_FAILURE;
    }

    destroy_mat(&result);
    
    return EXIT_SUCCESS;
}

int mode_partial(std::string table_filename, std::string tree_filename, 
                 std::string output_filename, std::string method_string,
                 bool vaw, double g_unifrac_alpha, bool bypass_tips, 
                 unsigned int nthreads, int start_stripe, int stop_stripe) {
    if(output_filename.empty()) {
        err("output filename missing");
        return EXIT_FAILURE;
    }

    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }

    if(tree_filename.empty()) {
        err("tree filename missing");
        return EXIT_FAILURE;
    }
    
    if(method_string.empty()) {
        err("method missing");
        return EXIT_FAILURE;
    }

    if(start_stripe < 0) {
        err("Starting stripe must be >= 0");
        return EXIT_FAILURE;
    }
    if(stop_stripe <= start_stripe) {
        err("In '--mode partial', the stop and start stripes must be specified, and the stop stripe must be > start stripe");
        return EXIT_FAILURE;
    }

    partial_mat_t *result = NULL;
    compute_status status;
    status = partial(table_filename.c_str(), tree_filename.c_str(), method_string.c_str(), 
                     vaw, g_unifrac_alpha, bypass_tips, nthreads, start_stripe, stop_stripe, &result);
    if(status != okay || result == NULL) {
        fprintf(stderr, "Compute failed in one_off with error code: %d\n", status);
        exit(EXIT_FAILURE);
    }
   
    io_status err = write_partial(output_filename.c_str(), result);
    destroy_partial_mat(&result);

    if(err != write_okay){
        fprintf(stderr, "Write failed: %s\n", err == open_error ? "could not open output" : "unknown error");
        return EXIT_FAILURE;
    } 

    return EXIT_SUCCESS;
}

int mode_one_off(std::string table_filename, std::string tree_filename, 
                 std::string output_filename, std::string method_string,
                 bool vaw, double g_unifrac_alpha, bool bypass_tips,
                 unsigned int nthreads) {
    if(output_filename.empty()) {
        err("output filename missing");
        return EXIT_FAILURE;
    }

    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }

    if(tree_filename.empty()) {
        err("tree filename missing");
        return EXIT_FAILURE;
    }
    
    if(method_string.empty()) {
        err("method missing");
        return EXIT_FAILURE;
    }

    mat_t *result = NULL;
    compute_status status;
    status = one_off(table_filename.c_str(), tree_filename.c_str(), method_string.c_str(), 
                     vaw, g_unifrac_alpha, bypass_tips, nthreads, &result);
    if(status != okay || result == NULL) {
        fprintf(stderr, "Compute failed in one_off with error code: %d\n", status);
        exit(EXIT_FAILURE);
    }
   
    write_mat(output_filename.c_str(), result);
    destroy_mat(&result);

    return EXIT_SUCCESS;
}

int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help")) {
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
    const std::string &partial_pattern = input.getCmdOption("--partial-pattern");
    const std::string &npartials = input.getCmdOption("--n-partials");

    if(nthreads_arg.empty()) {
        nthreads = 1;
    } else {
        nthreads = atoi(nthreads_arg.c_str());
    }
    
    bool vaw = input.cmdOptionExists("--vaw"); 
    bool bypass_tips = input.cmdOptionExists("-f");
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

    int n_partials;
    if(npartials.empty()) 
        n_partials = 1;
    else
        n_partials = atoi(npartials.c_str());

    if(mode_arg.empty() || mode_arg == "one-off")
        return mode_one_off(table_filename, tree_filename, output_filename, method_string, vaw, g_unifrac_alpha, bypass_tips, nthreads);
    else if(mode_arg == "partial")
        return mode_partial(table_filename, tree_filename, output_filename, method_string, vaw, g_unifrac_alpha, bypass_tips, nthreads, start_stripe, stop_stripe);
    else if(mode_arg == "merge-partial")
        return mode_merge_partial(output_filename, partial_pattern, nthreads);
    else if(mode_arg == "partial-report")
        return mode_partial_report(table_filename, n_partials);
    else 
        err("Unknown mode. Valid options are: one-off, partial, merge-partial");

    return EXIT_SUCCESS;
}


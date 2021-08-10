#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <glob.h>
#include <signal.h>
#include "api.hpp"
#include "cmd.hpp"
#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"

enum Format {format_invalid,format_ascii, format_hdf5_fp32, format_hdf5_fp64};

void usage() {
    std::cout << "usage: ssu -i <biom> -o <out.dm> -m [METHOD] -t <newick> [-n threads] [-a alpha] [-f]  [--vaw]" << std::endl;
    std::cout << "    [--mode [MODE]] [--start starting-stripe] [--stop stopping-stripe] [--partial-pattern <glob>]" << std::endl;
    std::cout << "    [--n-partials number_of_partitions] [--report-bare] [--format|-r out-mode]" << std::endl;
    std::cout << std::endl;
    std::cout << "    -i\t\tThe input BIOM table." << std::endl;
    std::cout << "    -t\t\tThe input phylogeny in newick." << std::endl;
    std::cout << "    -m\t\tThe method, [unweighted | weighted_normalized | weighted_unnormalized | generalized | unweighted_fp32 | weighted_normalized_fp32 | weighted_unnormalized_fp32 | generalized_fp32]." << std::endl;
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
    std::cout << "    \t\t    check-partial : Check partial UniFrac results." << std::endl;
    std::cout << "    --start\t[OPTIONAL] If mode==partial, the starting stripe." << std::endl;
    std::cout << "    --stop\t[OPTIONAL] If mode==partial, the stopping stripe." << std::endl;
    std::cout << "    --partial-pattern\t[OPTIONAL] If mode==merge-partial or check-partial, a glob pattern for partial outputs to merge." << std::endl;
    std::cout << "    --n-partials\t[OPTIONAL] If mode==partial-report, the number of partitions to compute." << std::endl;
    std::cout << "    --report-bare\t[OPTIONAL] If mode==partial-report, produce barebones output." << std::endl;
    std::cout << "    --format|-r\t[OPTIONAL]  Output format:" << std::endl;
    std::cout << "    \t\t    ascii : [DEFAULT] Original ASCII format." << std::endl;
    std::cout << "    \t\t    hfd5 : HFD5 format.  May be fp32 or fp64, depending on method." << std::endl;
    std::cout << "    \t\t    hdf5_fp32 : HFD5 format, using fp32 precision." << std::endl;
    std::cout << "    \t\t    hdf5_fp64 : HFD5 format, using fp64 precision." << std::endl;
    std::cout << "    --pcoa\t[OPTIONAL] Number of PCoA dimensions to compute (default: 10, do not compute if 0)" << std::endl;
    std::cout << "    --diskbuf\t[OPTIONAL] Use a disk buffer to reduce memory footprint. Provide path to a fast partition (ideally NVMe)." << std::endl;
    std::cout << std::endl;
    std::cout << "Citations: " << std::endl;
    std::cout << "    For UniFrac, please see:" << std::endl;
    std::cout << "        McDonald et al. Nature Methods 2018; DOI: 10.1038/s41592-018-0187-8" << std::endl;
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

const char* compute_status_messages[7] = {"No error.",
                                          "The tree file cannot be found.", 
                                          "The table file cannot be found.",
                                          "The table file contains an empty table.",
                                          "An unknown method was requested.", 
                                          "Table observation IDs are not a subset of the tree tips. This error can also be triggered if a node name contains a single quote (this is unlikely).",
                                          "Error creating the output."};


// https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
inline std::vector<std::string> glob(const std::string& pat){
    using namespace std;
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}


void err(std::string msg) {
    std::cerr << "ERROR: " << msg << std::endl << std::endl;
    usage();
}

int mode_partial_report(const std::string table_filename, unsigned int npartials, bool bare) {
    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }
    
    if(npartials < 1) {
        err("--n-partials cannot be < 1");
        exit(EXIT_FAILURE);
    }

    su::biom table = su::biom(table_filename.c_str());
    int total_stripes = (table.n_samples + 1) / 2;

    if(!bare) {
        std::cout << "Total samples: " << table.n_samples << std::endl;
        std::cout << "Total stripes: " << total_stripes << std::endl;
    }

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
            if(bare) 
                std::cout << start << "\t" << stop << std::endl;
            else
                std::cout << "Partition " << p << ", suggested start and stop: " << start << ", " << stop << std::endl;
            start = start + fullchunk;
        } else {
            stop = start + smallchunk;  // stripe end 
            if(bare) 
                std::cout << start << "\t" << stop << std::endl;
            else
                std::cout << "Partition " << p << ", suggested start and stop: " << start << ", " << stop << std::endl;
            start = start + smallchunk;
        }
    }

    return EXIT_SUCCESS;
} 

int mode_merge_partial_fp32(const char * output_filename, Format format_val, unsigned int pcoa_dims,
                            size_t partials_size, partial_dyn_mat_t* * partial_mats,
                            const char * mmap_dir) {
    mat_full_fp32_t *result = NULL;

    MergeStatus status = merge_partial_to_mmap_matrix_fp32(partial_mats, partials_size, mmap_dir, &result);

    if(status != merge_okay) {
        std::ostringstream msg;
        msg << "Unable to complete merge; err " << status;
        err(msg.str());
        return EXIT_FAILURE;
    }

    IOStatus iostatus;
    iostatus = write_mat_from_matrix_hdf5_fp32(output_filename, result, pcoa_dims);
    destroy_mat_full_fp32(&result);
    
    if(iostatus != write_okay) {
        std::ostringstream msg; 
        msg << "Unable to write; err " << iostatus;
        err(msg.str()); 
        return EXIT_FAILURE;
    }   
        
    return EXIT_SUCCESS;
}

int mode_merge_partial_fp64(const char * output_filename, Format format_val, unsigned int pcoa_dims,
                            size_t partials_size, partial_dyn_mat_t* * partial_mats,
                            const char * mmap_dir) {
    mat_full_fp64_t *result = NULL;

    MergeStatus status = merge_partial_to_mmap_matrix(partial_mats, partials_size, mmap_dir, &result);

    if(status != merge_okay) {
        std::ostringstream msg;
        msg << "Unable to complete merge; err " << status;
        err(msg.str());
        return EXIT_FAILURE;
    }

    IOStatus iostatus;
    if (format_val==format_hdf5_fp64) {
     iostatus = write_mat_from_matrix_hdf5(output_filename, result, pcoa_dims);
    } else {
     iostatus = write_mat_from_matrix(output_filename, result);
    }
    destroy_mat_full_fp64(&result);

    if(iostatus != write_okay) {
        std::ostringstream msg;
        msg << "Unable to write; err " << iostatus;
        err(msg.str());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int mode_merge_partial(const std::string &output_filename, Format format_val, unsigned int pcoa_dims,
                       const std::string &partial_pattern,
                       const std::string &mmap_dir) {
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
    
    std::vector<std::string> partials = glob(partial_pattern);
    partial_dyn_mat_t** partial_mats = (partial_dyn_mat_t**)malloc(sizeof(partial_dyn_mat_t*) * partials.size());
    for(size_t i = 0; i < partials.size(); i++) {
        IOStatus io_err = read_partial_header(partials[i].c_str(), &partial_mats[i]);
        if(io_err != read_okay) {
            std::ostringstream msg;
            msg << "Unable to parse file (" << partials[i] << "); err " << io_err;
            err(msg.str());
            return EXIT_FAILURE;
        }
    }

    const char * mmap_dir_c = mmap_dir.empty() ? NULL : mmap_dir.c_str();

    int status;
    if (format_val==format_hdf5_fp64) {
     status = mode_merge_partial_fp64(output_filename.c_str(), format_val, pcoa_dims, partials.size(), partial_mats, mmap_dir_c);
    } else if (format_val==format_hdf5_fp32) {
     status = mode_merge_partial_fp32(output_filename.c_str(), format_val, pcoa_dims, partials.size(), partial_mats, mmap_dir_c);
    } else {
     status = mode_merge_partial_fp64(output_filename.c_str(), format_val, pcoa_dims, partials.size(), partial_mats, mmap_dir_c);
    }

    for(size_t i = 0; i < partials.size(); i++) {
      destroy_partial_dyn_mat(&partial_mats[i]);
    }

    return status;
}

int mode_check_partial(const std::string &partial_pattern) {
    if(partial_pattern.empty()) {
        std::string msg("Partial file pattern missing. For instance, if your partial results\n" \
                        "are named 'ssu.unweighted.start0.partial', 'ssu.unweighted.start10.partial', \n" \
                        "etc, then a pattern of 'ssu.unweighted.start*.partial' would make sense");
        err(msg);
        return EXIT_FAILURE;
    }

    std::vector<std::string> partials = glob(partial_pattern);
    std::cout << "Processing " << partials.size() << " partial files." << std::endl;


    partial_dyn_mat_t** partial_mats = (partial_dyn_mat_t**)malloc(sizeof(partial_dyn_mat_t*) * partials.size());
    for(size_t i = 0; i < partials.size(); i++) {
        IOStatus io_err = read_partial_header(partials[i].c_str(), &partial_mats[i]);
        if(io_err != read_okay) {
            std::ostringstream msg;
            msg << "Unable to parse file (" << partials[i] << "); err " << io_err;
            err(msg.str());
            return EXIT_FAILURE;
        }
    }

    int status = validate_partial(partial_mats,partials.size());

    if (status==merge_okay) {
        std::cout << "Partials are ready to be merged" << std::endl;
        uint64_t n_samples_64 = partial_mats[0]->n_samples;
        uint64_t msize32 = sizeof(float) * n_samples_64 * n_samples_64;
        uint64_t msize64 = sizeof(double) * n_samples_64 * n_samples_64;

        std::cout << "n_samples = " << n_samples_64 << std::endl;
        std::cout << "Matrix RAM needs: fp32 = " << msize32/(1024.0*1024*1024) << "GB, fp64 = " << msize64/(1024.0*1024*1024) << "GB"<< std::endl;
    }

    for(size_t i = 0; i < partials.size(); i++) {
      destroy_partial_dyn_mat(&partial_mats[i]);
    }

    return status;
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
        fprintf(stderr, "Compute failed in partial: %s\n", compute_status_messages[status]);
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

int mode_one_off(const std::string &table_filename, const std::string &tree_filename, 
                 const std::string &output_filename, const std::string &format_str, Format format_val, 
                 const std::string &method_string, unsigned int pcoa_dims,
                 bool vaw, double g_unifrac_alpha, bool bypass_tips,
                 unsigned int nthreads, const std::string &mmap_dir) {
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

    compute_status status = okay;
    if (format_val==format_ascii) {
      mat_t *result = NULL;

      status = one_off(table_filename.c_str(), tree_filename.c_str(), method_string.c_str(), 
                       vaw, g_unifrac_alpha, bypass_tips, nthreads, &result);
      if(status != okay || result == NULL) {
        fprintf(stderr, "Compute failed in one_off: %s\n", compute_status_messages[status]);
        exit(EXIT_FAILURE);
      }
  
      IOStatus iostatus = write_mat(output_filename.c_str(), result);
      destroy_mat(&result);

      if(iostatus!=write_okay) {
        err("Failed to write output file.");
        status = output_error;
      }

    } else {
      const char * mmap_dir_c = mmap_dir.empty() ? NULL : mmap_dir.c_str();

      status = unifrac_to_file(table_filename.c_str(), tree_filename.c_str(), output_filename.c_str(),
                               method_string.c_str(), vaw, g_unifrac_alpha, bypass_tips, nthreads, format_str.c_str(),
                               pcoa_dims, mmap_dir_c);

      if (status != okay) {
        fprintf(stderr, "Compute failed in one_off: %s\n", compute_status_messages[status]);
      }
    }

    return (status==okay) ? EXIT_SUCCESS : EXIT_FAILURE;
}

void ssu_sig_handler(int signo) {
    if (signo == SIGUSR1) {
        printf("Status cannot be reported.\n");
    }
}

Format get_format(const std::string &format_string, const std::string &method_string) {
    Format format_val = format_invalid;
    if (format_string.empty()) {
        format_val = format_ascii;
    } else if (format_string == "ascii") {
        format_val = format_ascii;
    } else if (format_string == "hdf5_fp32") {
        format_val = format_hdf5_fp32;
    } else if (format_string == "hdf5_fp64") {
        format_val = format_hdf5_fp64;
    } else if (format_string == "hdf5") {
        if ((method_string=="unweighted_fp32") || (method_string=="weighted_normalized_fp32") || (method_string=="weighted_unnormalized_fp32") || (method_string=="generalized_fp32"))
           format_val = format_hdf5_fp32;
        else
           format_val = format_hdf5_fp64;
    }

    return format_val;
}

int main(int argc, char **argv){
    signal(SIGUSR1, ssu_sig_handler);
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help") || argc == 1) {
        usage();
        return EXIT_SUCCESS;
    }

    unsigned int nthreads;
    std::string table_filename = input.getCmdOption("-i");
    std::string tree_filename = input.getCmdOption("-t");
    std::string output_filename = input.getCmdOption("-o");
    std::string method_string = input.getCmdOption("-m");
    std::string nthreads_arg = input.getCmdOption("-n");
    std::string gunifrac_arg = input.getCmdOption("-a");
    std::string mode_arg = input.getCmdOption("--mode");
    std::string start_arg = input.getCmdOption("--start");
    std::string stop_arg = input.getCmdOption("--stop");
    std::string partial_pattern = input.getCmdOption("--partial-pattern");
    std::string npartials = input.getCmdOption("--n-partials");
    std::string report_bare = input.getCmdOption("--report-bare");
    std::string format_arg = input.getCmdOption("--format");
    std::string sformat_arg = input.getCmdOption("-r");
    std::string pcoa_arg = input.getCmdOption("--pcoa");
    std::string diskbuf_arg = input.getCmdOption("--diskbuf");

    if(nthreads_arg.empty()) {
        nthreads = 1;
    } else {
        nthreads = atoi(nthreads_arg.c_str());
    }
    
    bool vaw = input.cmdOptionExists("--vaw"); 
    bool bare = input.cmdOptionExists("--report-bare"); 
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

    if(n_partials<1) {
        err("--n-partials cannot be < 1");
        return EXIT_FAILURE;
    }
     if(n_partials>1000000000) {
        err("--n-partials cannot be > 1G");
        return EXIT_FAILURE;
    }

    Format format_val = format_invalid;
    if(!format_arg.empty()) {
      format_val = get_format(format_arg,method_string);
    } else {
      format_val = get_format(sformat_arg,method_string);
      format_arg=sformat_arg; // easier to use a single variable
    }
    if(format_val==format_invalid) {
        err("Invalid format, must be one of ascii|hdf5|hdf5_fp32|hdf5_fp64");
        return EXIT_FAILURE;
    }

    unsigned int pcoa_dims;
    if(pcoa_arg.empty())
        pcoa_dims = 10;
    else
        pcoa_dims = atoi(pcoa_arg.c_str());


    if(mode_arg.empty() || mode_arg == "one-off")
        return mode_one_off(table_filename, tree_filename, output_filename, format_arg, format_val, method_string, pcoa_dims, vaw, g_unifrac_alpha, bypass_tips, nthreads, diskbuf_arg);
    else if(mode_arg == "partial")
        return mode_partial(table_filename, tree_filename, output_filename, method_string, vaw, g_unifrac_alpha, bypass_tips, nthreads, start_stripe, stop_stripe);
    else if(mode_arg == "merge-partial")
        return mode_merge_partial(output_filename, format_val, pcoa_dims, partial_pattern, diskbuf_arg);
    else if(mode_arg == "check-partial")
        return mode_check_partial(partial_pattern);
    else if(mode_arg == "partial-report")
        return mode_partial_report(table_filename, uint32_t(n_partials), bare);
    else 
        err("Unknown mode. Valid options are: one-off, partial, merge-partial, check-partial, partial-report");

    return EXIT_SUCCESS;
}


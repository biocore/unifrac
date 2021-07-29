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
    std::cout << "usage: faithpd -i <biom> -t <newick> -o <out.txt>" << std::endl;
    std::cout << std::endl;
    std::cout << "    -i\t\tThe input BIOM table." << std::endl;
    std::cout << "    -t\t\tThe input phylogeny in newick." << std::endl;
    std::cout << "    -o\t\tThe output series." << std::endl;
    std::cout << std::endl;
    std::cout << "Citations: " << std::endl;
    std::cout << "    For Faith's PD, please see:" << std::endl;
    std::cout << "        Faith Biological Conservation 1992; DOI: 10.1016/0006-3207(92)91201-3" << std::endl;
    std::cout << std::endl;

}

const char* compute_status_messages[7] = {"No error.",
                                          "The tree file cannot be found.",
                                          "The table file cannot be found.",
                                          "The table file contains an empty table.",
                                          "An unknown method was requested.",
                                          "Table observation IDs are not a subset of the tree tips. This error can also be triggered if a node name contains a single quote (this is unlikely).",
                                          "Error creating the output."};

void err(std::string msg) {
    std::cerr << "ERROR: " << msg << std::endl << std::endl;
    usage();
}

int faith_cli_one_off(std::string table_filename, std::string tree_filename,
                 std::string output_filename) {
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

    r_vec *result = NULL;
    compute_status status;
    status = faith_pd_one_off(table_filename.c_str(), tree_filename.c_str(), &result);
    if(status != okay || result == NULL) {
        fprintf(stderr, "Compute failed in faith_pd_one_off: %s\n", compute_status_messages[status]);
        exit(EXIT_FAILURE);
    }

    write_vec(output_filename.c_str(), result);
    destroy_results_vec(&result);

    return EXIT_SUCCESS;
}

int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help") || argc == 1) {
        usage();
        return EXIT_SUCCESS;
    }

    const std::string &table_filename = input.getCmdOption("-i");
    const std::string &tree_filename = input.getCmdOption("-t");
    const std::string &output_filename = input.getCmdOption("-o");

    faith_cli_one_off(table_filename, tree_filename, output_filename);

    return EXIT_SUCCESS;
}

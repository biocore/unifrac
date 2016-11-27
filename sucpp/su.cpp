#include <iostream>
#include <fstream>
#include <string>
#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include <iomanip>

int main(int argc, char** argv) {
    std::ifstream ifs(argv[1]);
    std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    
    su::BPTree tree = su::BPTree(content);
    su::biom table = su::biom(argv[2]);
    su::Method method;

    if(strcmp(argv[3], "unweighted") == 0)
        method = su::unweighted;
    else if(strcmp(argv[3], "weighted_normalized") == 0) 
        method = su::weighted_normalized;
    else if(strcmp(argv[3], "weighted_unnormalized") == 0)
        method = su::weighted_unnormalized;
    else {
        std::cerr << "Unknown method: " << argv[3] << std::endl;
        return EXIT_FAILURE;
    }

    double **dm = unifrac(table, tree, method);
    for(unsigned int i = 0; i < table.n_samples; i++)
        std::cout << "\t" << table.sample_ids[i];
    std::cout << std::endl;
    for(unsigned int i = 0; i < table.n_samples; i++) {
        std::cout << table.sample_ids[i];
        for(unsigned int j = 0; j < table.n_samples; j++)
            std::cout << std::setprecision(16) << "\t" << dm[i][j];
        std::cout << std::endl;
    }
    return EXIT_SUCCESS;
}

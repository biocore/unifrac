#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace su {

    class BPTree {
        public:
            std::vector<double> lengths;
            std::vector<std::string> names;
            uint32_t nparens;

            BPTree(std::string newick);
            ~BPTree();

            uint32_t postorderselect(uint32_t i);
            bool isleaf(uint32_t i);
            uint32_t leftchild(uint32_t i);
            uint32_t rightchild(uint32_t i);
            uint32_t rightsibling(uint32_t i);
            std::vector<bool> get_structure();
            std::vector<uint32_t> get_openclose();

            void print() {
                for(auto c = structure.begin(); c != structure.end(); c++) {
                    if(*c)
                        std::cout << "1";
                    else
                        std::cout << "0";
                }
                std::cout << std::endl;
            }
        private:
            std::vector<bool> structure;
            std::vector<uint32_t> openclose;
            std::vector<uint32_t> select_0_index;
            std::vector<uint32_t> select_1_index;

            void index_select();
            void newick_to_bp(std::string newick);
            void newick_to_metadata(std::string newick);
            void structure_to_openclose();
            void set_node_metadata(unsigned int open_idx, std::string &token);
            bool is_structure_character(char c);
            inline uint32_t open(uint32_t i);
            inline uint32_t close(uint32_t i);
            std::string tokenize(std::string::iterator &start, const std::string::iterator &end);
    };
}

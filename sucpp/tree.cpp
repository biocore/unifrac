#include "tree.hpp"
#include <stack>

using namespace su;

BPTree::BPTree(std::string newick) {
    openclose = std::vector<uint32_t>();
    lengths = std::vector<double>();
    names = std::vector<std::string>();
    select_0_index = std::vector<uint32_t>();
    select_1_index = std::vector<uint32_t>();
    structure = std::vector<bool>();
    structure.reserve(500000);

    // three pass for parse. not ideal, but easier to map from IOW code    
    newick_to_bp(newick);

    // resize is correct here as we are not performing a push_back
    openclose.resize(nparens);
    lengths.resize(nparens);
    names.resize(nparens);
    select_0_index.resize(nparens / 2);
    select_1_index.resize(nparens / 2);

    structure_to_openclose();
    newick_to_metadata(newick);
    index_select();
}

BPTree::~BPTree() {
    //bit_array_free(B);
}

void BPTree::index_select() {
    unsigned int idx = 0;
    auto i = structure.begin();
    auto k0 = select_0_index.begin();
    auto k1 = select_1_index.begin();

    for(; i != structure.end(); i++, idx++ ) {
        if(*i) 
            *(k1++) = idx;
        else
            *(k0++) = idx;
    }
}

uint32_t BPTree::postorderselect(uint32_t k) { 
    return open(select_0_index[k]);
}

inline uint32_t BPTree::open(uint32_t i) {
    return structure[i] ? i : openclose[i];
}

inline uint32_t BPTree::close(uint32_t i) {
    return structure[i] ? openclose[i] : i;
}

bool BPTree::isleaf(unsigned int idx) {
    return (structure[idx] && !structure[idx + 1]);
}

uint32_t BPTree::leftchild(uint32_t i) {
    // aka fchild
    if(isleaf(i))
        return 0;  // this is awkward, using 0 which is root, but a root cannot be a child. edge case
    else
        return i + 1;
}

uint32_t BPTree::rightchild(uint32_t i) {
    // aka lchild
    if(isleaf(i))
        return 0;  // this is awkward, using 0 which is root, but a root cannot be a child. edge case
    else
        return open(close(i) - 1);
}

uint32_t BPTree::rightsibling(uint32_t i) {
    // aka nsibling
    uint32_t position = close(i) + 1;
    if(position >= nparens)
        return 0;  // will return 0 if no sibling as root cannot have a sibling
    else if(structure[position])
        return position;
    else 
        return 0;
}

void BPTree::newick_to_bp(std::string newick) {
    char last_structure;
    bool potential_single_descendent = false;
    int count = 0;
    for(auto c = newick.begin(); c != newick.end(); c++) {
        switch(*c) {
            case '(':
                // opening of a node
                count++;
                structure.push_back(true);
                last_structure = *c;
                potential_single_descendent = true;
                break;
            case ')':
                // closing of a node
                if(potential_single_descendent || (last_structure == ',')) {
                    // we have a single descendent or a last child (i.e. ",)" scenario)
                    count += 3;
                    structure.push_back(true);
                    structure.push_back(false);
                    structure.push_back(false);
                    potential_single_descendent = false;
                } else {
                    // it is possible still to have a single descendent in the case of 
                    // multiple single descendents (e.g., (...()...) )
                    count += 1;
                    structure.push_back(false);
                }
                last_structure = *c;
                break;
            case ',':
                if(last_structure != ')') {
                    // we have a new tip
                    count += 2;
                    structure.push_back(true);
                    structure.push_back(false);
                }
                potential_single_descendent = false;
                last_structure = *c;
                break;
            default:
                break;
        }
    }
    nparens = structure.size();
}


void BPTree::structure_to_openclose() {
    std::stack<unsigned int> oc;
    unsigned int open_idx;
    unsigned int i = 0;

    for(auto it = structure.begin(); it != structure.end(); it++, i++) {
        if(*it) {
            oc.push(i);
        } else {
            open_idx = oc.top();
            oc.pop();
            openclose[i] = open_idx;
            openclose[open_idx] = i;
        }
    }
}


void BPTree::newick_to_metadata(std::string newick) {
    std::string::iterator start = newick.begin();
    std::string::iterator end = newick.end();
    std::string token;
    char last_structure = '\0';

    unsigned int structure_idx = 0;
    unsigned int lag = 0;
    unsigned int open_idx;

    while(start != end) {
        token = tokenize(start, end);
        // this sucks. 
        if(token.length() == 1 && is_structure_character(token[0])) {
            switch(token[0]) {
                case '(':
                    structure_idx++;
                    break;
                case ')':
                case ',':
                    structure_idx++;
                    if(last_structure == ')')
                        lag++;
                    break;
            }
        } else {
            // puts us on the corresponding closing parenthesis
            structure_idx += lag;
            lag = 0;

            open_idx = open(structure_idx);
            set_node_metadata(open_idx, token);
            // std::cout << structure_idx << " <-> " << open_idx << " " << token << std::endl;
            // make sure to advance an extra position if we are a leaf as the
            // as a leaf is by definition a 10, and doing a single advancement
            // would put the structure to token mapping out of sync
            if(isleaf(open_idx))
                structure_idx += 2;
            else
                structure_idx += 1;
            
        }
        last_structure = token[0];  
    }
}

void BPTree::set_node_metadata(unsigned int open_idx, std::string &token) {
    double length = 0.0;
    std::string name = std::string();
    unsigned int colon_idx = token.find_last_of(':');
    
    if(colon_idx == 0)
        length = std::stof(token.substr(1));
    else if(colon_idx < token.length()) {
        name = token.substr(0, colon_idx);
        length = std::stof(token.substr(colon_idx + 1));
    } else 
        name = token;
    
    names[open_idx] = name;
    lengths[open_idx] = length;
}

inline bool BPTree::is_structure_character(char c) {
    return (c == '(' || c == ')' || c == ',' || c == ';');
}

std::string BPTree::tokenize(std::string::iterator &start, const std::string::iterator &end) {
    bool inquote = false;
    bool isquote = false;
    char c;
    std::string token;
    
    do {
        c = *start;
        start++;

        isquote = (c == '"' || c == '\'');
        
        if(inquote && isquote) {
            inquote = false;
            continue;
        } else if(!inquote && isquote) {
            inquote = true;
            continue;
        }
    
        if(is_structure_character(c) && !inquote) {
            if(token.length() == 0)
                token.push_back(c);
            break;
        }

        token.push_back(c);
            

    } while(start != end);
    
    return token;
}   

std::vector<bool> BPTree::get_structure() {
    return structure;
}

std::vector<uint32_t> BPTree::get_openclose() {
    return openclose;
}

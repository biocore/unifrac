#include "tree.hpp"
#include <stack>
#include <algorithm>

using namespace su;

BPTree::BPTree(std::string newick) {
    openclose = std::vector<uint32_t>();
    lengths = std::vector<double>();
    names = std::vector<std::string>();
    excess = std::vector<uint32_t>();

    select_0_index = std::vector<uint32_t>();
    select_1_index = std::vector<uint32_t>();
    structure = std::vector<bool>();
    structure.reserve(500000);  // a fair sized tree... avoid reallocs, and its not _that_ much waste if this is wrong

    // three pass for parse. not ideal, but easier to map from IOW code    
    newick_to_bp(newick);

    // resize is correct here as we are not performing a push_back
    openclose.resize(nparens);
    lengths.resize(nparens);
    names.resize(nparens);
    select_0_index.resize(nparens / 2);
    select_1_index.resize(nparens / 2);
    excess.resize(nparens);

    structure_to_openclose();
    newick_to_metadata(newick);
    index_and_cache();
}

BPTree::BPTree(std::vector<bool> input_structure, std::vector<double> input_lengths, std::vector<std::string> input_names) {
    structure = input_structure;
    lengths = input_lengths;
    names = input_names;
    
    nparens = structure.size();

    openclose = std::vector<uint32_t>();
    select_0_index = std::vector<uint32_t>();
    select_1_index = std::vector<uint32_t>();
    openclose.resize(nparens);
    select_0_index.resize(nparens / 2);
    select_1_index.resize(nparens / 2);
    excess.resize(nparens);

    structure_to_openclose();
    index_and_cache();
}

BPTree BPTree::mask(std::vector<bool> topology_mask, std::vector<double> in_lengths) {
    
    std::vector<bool> new_structure = std::vector<bool>();
    std::vector<double> new_lengths = std::vector<double>();
    std::vector<std::string> new_names = std::vector<std::string>();

    uint32_t count = 0;
    for(auto i = topology_mask.begin(); i != topology_mask.end(); i++) {
        if(*i)
            count++;
    }

    new_structure.resize(count);
    new_lengths.resize(count);
    new_names.resize(count);

    auto mask_it = topology_mask.begin();
    auto base_it = this->structure.begin();
    uint32_t new_idx = 0;
    uint32_t old_idx = 0;
    for(; mask_it != topology_mask.end(); mask_it++, base_it++, old_idx++) {
        if(*mask_it) {
            new_structure[new_idx] = this->structure[old_idx];
            new_lengths[new_idx] = in_lengths[old_idx];
            new_names[new_idx] = this->names[old_idx];
            new_idx++;
        }
    }
    
    return BPTree(new_structure, new_lengths, new_names);
}

std::unordered_set<std::string> BPTree::get_tip_names() {
    std::unordered_set<std::string> observed;
	
    for(unsigned int i = 0; i < this->nparens; i++) {
        if(this->isleaf(i)) {
            observed.insert(this->names[i]);
        }
    }

    return observed;
}

BPTree BPTree::shear(std::unordered_set<std::string> to_keep) {
    std::vector<bool> shearmask = std::vector<bool>(this->nparens);
    int32_t p;

	for(unsigned int i = 0; i < this->nparens; i++) {
        if(this->isleaf(i) && to_keep.count(this->names[i]) > 0) {
            shearmask[i] = true;
            shearmask[i+1] = true;

            p = this->parent(i);
            while(p != -1 && !shearmask[p]) {
                shearmask[p] = true;
                shearmask[this->close(p)] = true;
                p = this->parent(p);
            }
        }
    }
    return this->mask(shearmask, this->lengths);
}

BPTree BPTree::collapse() {
    std::vector<bool> collapsemask = std::vector<bool>(this->nparens);
    std::vector<double> new_lengths = std::vector<double>(this->lengths);

    uint32_t current, first, last;

    for(uint32_t i = 0; i < this->nparens / 2; i++) {
        current = this->preorderselect(i);

        if(this->isleaf(current) or (current == 0)) {  // 0 == root
            collapsemask[current] = true;
            collapsemask[this->close(current)] = true;
        } else {
            first = this->leftchild(current);
            last = this->rightchild(current);

            if(first == last) {
                new_lengths[first] = new_lengths[first] + new_lengths[current];
            } else {
                collapsemask[current] = true;
                collapsemask[this->close(current)] = true;
            }
        }
    }

    return this->mask(collapsemask, new_lengths);
}
   /*
        mask = bit_array_create(self.B.size)
        bit_array_set_bit(mask, self.root())
        bit_array_set_bit(mask, self.close(self.root()))

        new_lengths = self._lengths.copy()
        new_lengths_ptr = <DOUBLE_t*>new_lengths.data

        with nogil:
            for i in range(n):
                current = self.preorderselect(i)

                if self.isleaf(current):
                    bit_array_set_bit(mask, current)
                    bit_array_set_bit(mask, self.close(current))
                else:
                    first = self.fchild(current)
                    last = self.lchild(current)

                    if first == last:
                        new_lengths_ptr[first] = new_lengths_ptr[first] + \
                                new_lengths_ptr[current]
                    else:
                        bit_array_set_bit(mask, current)
                        bit_array_set_bit(mask, self.close(current))

        new_bp = self._mask_from_self(mask, new_lengths)
        bit_array_free(mask)
        return new_bp
*/ 


BPTree::~BPTree() {
}

void BPTree::index_and_cache() {
    // should probably do the open/close in here too
    unsigned int idx = 0;
    auto i = structure.begin();
    auto k0 = select_0_index.begin();
    auto k1 = select_1_index.begin();
    auto e_it = excess.begin();
    unsigned int e = 0;  

    for(; i != structure.end(); i++, idx++ ) {
        if(*i) {
            *(k1++) = idx;
            *(e_it++) = ++e;
        }
        else {
            *(k0++) = idx;
            *(e_it++) = --e;
        }
    }
}

uint32_t BPTree::postorderselect(uint32_t k) const { 
    return open(select_0_index[k]);
}

uint32_t BPTree::preorderselect(uint32_t k) const {
    return select_1_index[k];
}

inline uint32_t BPTree::open(uint32_t i) const {
    return structure[i] ? i : openclose[i];
}

inline uint32_t BPTree::close(uint32_t i) const {
    return structure[i] ? openclose[i] : i;
}

bool BPTree::isleaf(unsigned int idx) const {
    return (structure[idx] && !structure[idx + 1]);
}

uint32_t BPTree::leftchild(uint32_t i) const {
    // aka fchild
    if(isleaf(i))
        return 0;  // this is awkward, using 0 which is root, but a root cannot be a child. edge case
    else
        return i + 1;
}

uint32_t BPTree::rightchild(uint32_t i) const {
    // aka lchild
    if(isleaf(i))
        return 0;  // this is awkward, using 0 which is root, but a root cannot be a child. edge case
    else
        return open(close(i) - 1);
}

uint32_t BPTree::rightsibling(uint32_t i) const {
    // aka nsibling
    uint32_t position = close(i) + 1;
    if(position >= nparens)
        return 0;  // will return 0 if no sibling as root cannot have a sibling
    else if(structure[position])
        return position;
    else 
        return 0;
}

int32_t BPTree::parent(uint32_t i) const {
    return enclose(i);
}

int32_t BPTree::enclose(uint32_t i) const {
    if(structure[i])
        return bwd(i, -2) + 1;
    else
        return bwd(i - 1, -2) + 1; 
}

int32_t BPTree::bwd(uint32_t i, int d) const {
    uint32_t target_excess = excess[i] + d;
    for(int current_idx = i - 1; current_idx >= 0; current_idx--) {
        if(excess[current_idx] == target_excess)
            return current_idx;
    }
    return -1;
}

void BPTree::newick_to_bp(std::string newick) {
    char last_structure;
    bool potential_single_descendent = false;
    int count = 0;
    bool in_quote = false;
    for(auto c = newick.begin(); c != newick.end(); c++) {
        if(*c == '\'') 
            in_quote = !in_quote;

        if(in_quote)
            continue;

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
// trim from end
// from http://stackoverflow.com/a/217605
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}


//// WEIRDNESS. THIS SOLVES IT WITH THE RTRIM. ISOLATE, MOVE TO CONSTRUCTOR.
void BPTree::newick_to_metadata(std::string newick) {
    newick = rtrim(newick);
    
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

inline bool BPTree::is_structure_character(char c) const {
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
        
        if(c == '\n') {
            continue;
        }

        isquote = c == '\'';
        
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


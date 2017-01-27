#include <iostream>
#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include <cmath>
#include <unordered_set>
#include <string.h>

/*
 * test harness adapted from 
 * https://github.com/noporpoise/BitArray/blob/master/dev/bit_array_test.c
 */
const char *suite_name;
char suite_pass;
int suites_run = 0, suites_failed = 0, suites_empty = 0;
int tests_in_suite = 0, tests_run = 0, tests_failed = 0;

#define QUOTE(str) #str
#define ASSERT(x) {tests_run++; tests_in_suite++; if(!(x)) \
    { fprintf(stderr, "failed assert [%s:%i] %s\n", __FILE__, __LINE__, QUOTE(x)); \
      suite_pass = 0; tests_failed++; }}

void SUITE_START(const char *name) {
  suite_pass = 1;
  suite_name = name;
  suites_run++;
  tests_in_suite = 0;
}

void SUITE_END() {
  printf("Testing %s ", suite_name);
  size_t suite_i;
  for(suite_i = strlen(suite_name); suite_i < 80-8-5; suite_i++) printf(".");
  printf("%s\n", suite_pass ? " pass" : " fail");
  if(!suite_pass) suites_failed++;
  if(!tests_in_suite) suites_empty++;
}
/*
 *  End adapted code
 */

std::vector<bool> _bool_array_to_vector(bool *arr, unsigned int n) {
    std::vector<bool> vec;

    for(unsigned int i = 0; i < n; i++)
        vec.push_back(arr[i]);

    return vec;
}

std::vector<uint32_t> _uint32_array_to_vector(uint32_t *arr, unsigned int n) {
    std::vector<uint32_t> vec;

    for(unsigned int i = 0; i < n; i++)
        vec.push_back(arr[i]);

    return vec;
}

std::vector<double> _double_array_to_vector(double *arr, unsigned int n) {
    std::vector<double> vec;

    for(unsigned int i = 0; i < n; i++)
        vec.push_back(arr[i]);

    return vec;
}

std::vector<std::string> _string_array_to_vector(std::string *arr, unsigned int n) {
    std::vector<std::string> vec;

    for(unsigned int i = 0; i < n; i++)
        vec.push_back(arr[i]);

    return vec;
}

bool vec_almost_equal(std::vector<double> a, std::vector<double> b) {
    if(a.size() != b.size()) {
        return false;
    }
    for(unsigned int i = 0; i < a.size(); i++) {
        if(!(fabs(a[i] - b[i]) < 0.000001)) {  // sufficient given the tests
            return false;
        }
    }
    return true;
}


void test_bptree_constructor_simple() {
    SUITE_START("bptree constructor simple");
                                //01234567
                                //11101000
    su::BPTree tree = su::BPTree("(('123:foo; bar':1,b:2)c);");

    unsigned int exp_nparens = 8;

    std::vector<bool> exp_structure;
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(false);
    exp_structure.push_back(false);

    std::vector<uint32_t> exp_openclose;
    exp_openclose.push_back(7);
    exp_openclose.push_back(6);
    exp_openclose.push_back(3);
    exp_openclose.push_back(2);
    exp_openclose.push_back(5);
    exp_openclose.push_back(4);
    exp_openclose.push_back(1);
    exp_openclose.push_back(0);

    std::vector<std::string> exp_names;
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("c"));
    exp_names.push_back(std::string("123:foo; bar"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("b"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());

    std::vector<double> exp_lengths;
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(1.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(2.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(tree.get_openclose() == exp_openclose);
    ASSERT(tree.lengths == exp_lengths);
    ASSERT(tree.names == exp_names);

    SUITE_END();
}

void test_bptree_constructor_from_existing() {
    SUITE_START("bptree constructor from_existing");
                                //01234567
                                //11101000
    su::BPTree existing = su::BPTree("(('123:foo; bar':1,b:2)c);");
    su::BPTree tree = su::BPTree(existing.get_structure(), existing.lengths, existing.names);

    unsigned int exp_nparens = 8;

    std::vector<bool> exp_structure;
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(false);
    exp_structure.push_back(false);

    std::vector<uint32_t> exp_openclose;
    exp_openclose.push_back(7);
    exp_openclose.push_back(6);
    exp_openclose.push_back(3);
    exp_openclose.push_back(2);
    exp_openclose.push_back(5);
    exp_openclose.push_back(4);
    exp_openclose.push_back(1);
    exp_openclose.push_back(0);

    std::vector<std::string> exp_names;
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("c"));
    exp_names.push_back(std::string("123:foo; bar"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("b"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());

    std::vector<double> exp_lengths;
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(1.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(2.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(tree.get_openclose() == exp_openclose);
    ASSERT(tree.lengths == exp_lengths);
    ASSERT(tree.names == exp_names);

    SUITE_END();
}

void test_bptree_mask() {
    SUITE_START("bptree mask");
                                //01234567
                                //11101000
                                //111000
    std::vector<bool> mask = {true, true, true, true, false, false, true, true};
    su::BPTree base = su::BPTree("(('123:foo; bar':1,b:2)c);");
    su::BPTree tree = base.mask(mask, base.lengths);
    unsigned int exp_nparens = 6;

    std::vector<bool> exp_structure;
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(false);
    exp_structure.push_back(false);

    std::vector<uint32_t> exp_openclose;
    exp_openclose.push_back(5);
    exp_openclose.push_back(4);
    exp_openclose.push_back(3);
    exp_openclose.push_back(2);
    exp_openclose.push_back(1);
    exp_openclose.push_back(0);

    std::vector<std::string> exp_names;
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("c"));
    exp_names.push_back(std::string("123:foo; bar"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());

    std::vector<double> exp_lengths;
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(1.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(tree.get_openclose() == exp_openclose);
    ASSERT(tree.lengths == exp_lengths);
    ASSERT(tree.names == exp_names);

    SUITE_END();
}

void test_bptree_constructor_single_descendent() {
    SUITE_START("bptree constructor single descendent");

    su::BPTree tree = su::BPTree("(((a)b)c,((d)e)f,g)r;");
    
    unsigned int exp_nparens = 16;

    bool structure_arr[] = {1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0}; 
    std::vector<bool> exp_structure = _bool_array_to_vector(structure_arr, exp_nparens);

    double length_arr[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> exp_lengths = _double_array_to_vector(length_arr, exp_nparens);

    std::string names_arr[] = {"r", "c", "b", "a", "", "", "", "f", "e", "d", "", "", "", "g", "", ""};
    std::vector<std::string> exp_names = _string_array_to_vector(names_arr, exp_nparens);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(vec_almost_equal(tree.lengths, exp_lengths));
    ASSERT(tree.names == exp_names);
    
    SUITE_END();
}

void test_bptree_constructor_complex() {
    SUITE_START("bp tree constructor complex");
    su::BPTree tree = su::BPTree("(((a:1,b:2.5)c:6,d:8,(e),(f,g,(h:1,i:2)j:1)k:1.2)l,m:2)r;");
    
    unsigned int exp_nparens = 30;

    bool structure_arr[] = {1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0};
    std::vector<bool> exp_structure = _bool_array_to_vector(structure_arr, exp_nparens);

    double length_arr[] = {0, 0, 6, 1, 0, 2.5, 0, 0, 8, 0, 0, 0, 0, 0, 1.2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 2, 0, 0};
    std::vector<double> exp_lengths = _double_array_to_vector(length_arr, exp_nparens);

    std::string names_arr[] = {"r", "l", "c", "a", "", "b", "", "", "d", "", "", "e", "", "", "k", "f", "", "g", "", "j", "h", "", "i", "", "", "", "", "m", "", ""};
    std::vector<std::string> exp_names = _string_array_to_vector(names_arr, exp_nparens);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(vec_almost_equal(tree.lengths, exp_lengths));
    ASSERT(tree.names == exp_names);
    SUITE_END();
}

void test_bptree_constructor_semicolon() {
    SUITE_START("bp tree constructor semicolon");
    su::BPTree tree = su::BPTree("((a,(b,c):5)'d','e; foo':10,((f))g)r;");
    
    unsigned int exp_nparens = 20;

    bool structure_arr[] = {1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0};
    std::vector<bool> exp_structure = _bool_array_to_vector(structure_arr, exp_nparens);

    double length_arr[] = {0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> exp_lengths = _double_array_to_vector(length_arr, exp_nparens);

    std::string names_arr[] = {"r", "d", "a", "", "", "b", "", "c", "", "", "", "e; foo", "", "g", "", "f", "", "", "", ""};
    std::vector<std::string> exp_names = _string_array_to_vector(names_arr, exp_nparens);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(vec_almost_equal(tree.lengths, exp_lengths));
    ASSERT(tree.names == exp_names);
    SUITE_END();
}

void test_bptree_constructor_edgecases() {
    SUITE_START("bp tree constructor edgecases");

    su::BPTree tree1 = su::BPTree("((a,b));");
    bool structure_arr1[] = {1, 1, 1, 0, 1, 0, 0, 0};
    std::vector<bool> exp_structure1 = _bool_array_to_vector(structure_arr1, 8);

    su::BPTree tree2 = su::BPTree("(a);");
    bool structure_arr2[] = {1, 1, 0, 0};
    std::vector<bool> exp_structure2 = _bool_array_to_vector(structure_arr2, 4);

    su::BPTree tree3 = su::BPTree("();");
    bool structure_arr3[] = {1, 1, 0, 0};
    std::vector<bool> exp_structure3 = _bool_array_to_vector(structure_arr3, 4);

    su::BPTree tree4 = su::BPTree("((a,b),c);");
    bool structure_arr4[] = {1, 1, 1, 0, 1, 0, 0, 1, 0, 0};
    std::vector<bool> exp_structure4 = _bool_array_to_vector(structure_arr4, 10);

    su::BPTree tree5 = su::BPTree("(a,(b,c));");
    bool structure_arr5[] = {1, 1, 0, 1, 1, 0, 1, 0, 0, 0};
    std::vector<bool> exp_structure5 = _bool_array_to_vector(structure_arr5, 10);

    ASSERT(tree1.get_structure() == exp_structure1);
    ASSERT(tree2.get_structure() == exp_structure2);
    ASSERT(tree3.get_structure() == exp_structure3);
    ASSERT(tree4.get_structure() == exp_structure4);
    ASSERT(tree5.get_structure() == exp_structure5);

    SUITE_END();
}

void test_bptree_postorder() {
    SUITE_START("postorderselect");
    
    // fig1 from https://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    su::BPTree tree = su::BPTree("((3,4,(6)5)2,7,((10,100)9)8)1;");
    uint32_t exp[] = {2, 4, 7, 6, 1, 11, 15, 17, 14, 13, 0};
    uint32_t obs[tree.nparens / 2];
    
    for(int i = 0; i < (tree.nparens / 2); i++)
        obs[i] = tree.postorderselect(i);

    std::vector<uint32_t> exp_v = _uint32_array_to_vector(exp, tree.nparens / 2);
    std::vector<uint32_t> obs_v = _uint32_array_to_vector(obs, tree.nparens / 2);
    
    ASSERT(obs_v == exp_v);
    SUITE_END();
}

void test_bptree_preorder() {
    SUITE_START("preorderselect");
    
    // fig1 from https://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    su::BPTree tree = su::BPTree("((3,4,(6)5)2,7,((10,100)9)8)1;");
    uint32_t exp[] = {0, 1, 2, 4, 6, 7, 11, 13, 14, 15, 17};
    uint32_t obs[tree.nparens / 2];
    
    for(int i = 0; i < (tree.nparens / 2); i++)
        obs[i] = tree.preorderselect(i);

    std::vector<uint32_t> exp_v = _uint32_array_to_vector(exp, tree.nparens / 2);
    std::vector<uint32_t> obs_v = _uint32_array_to_vector(obs, tree.nparens / 2);
    
    ASSERT(obs_v == exp_v);
    SUITE_END();
}

void test_bptree_parent() {
    SUITE_START("parent");
    
    // fig1 from https://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    su::BPTree tree = su::BPTree("((3,4,(6)5)2,7,((10,100)9)8)1;");
    uint32_t exp[] = {0, 1, 1, 1, 1, 1, 6, 6, 1, 0, 0, 0, 0, 13, 14, 14, 14, 14, 13, 0};
    
    // all the -2 and +1 garbage is to avoid testing the root.
    uint32_t obs[tree.nparens - 2];
    
    for(int i = 0; i < (tree.nparens) - 2; i++)
        obs[i] = tree.parent(i+1);

    std::vector<uint32_t> exp_v = _uint32_array_to_vector(exp, tree.nparens - 2);
    std::vector<uint32_t> obs_v = _uint32_array_to_vector(obs, tree.nparens - 2);
   
    ASSERT(obs_v == exp_v);
    SUITE_END();
}

void test_biom_constructor() {
    SUITE_START("biom constructor");

    su::biom table = su::biom("test.biom");
    uint32_t exp_n_samples = 6;
    uint32_t exp_n_obs = 5;

    std::string sids[] = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};
    std::vector<std::string> exp_sids = _string_array_to_vector(sids, exp_n_samples);
    
    std::string oids[] = {"GG_OTU_1", "GG_OTU_2","GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
    std::vector<std::string> exp_oids = _string_array_to_vector(oids, exp_n_obs);

    uint32_t s_indptr[] = {0, 2, 5, 9, 11, 12, 15};
    std::vector<uint32_t> exp_s_indptr = _uint32_array_to_vector(s_indptr, exp_n_samples + 1);
    
    uint32_t o_indptr[] = {0, 1, 6, 9, 13, 15};
    std::vector<uint32_t> exp_o_indptr = _uint32_array_to_vector(o_indptr, exp_n_obs + 1);
  
    uint32_t exp_nnz = 15;

    ASSERT(table.n_samples == exp_n_samples);
    ASSERT(table.n_obs == exp_n_obs);
    ASSERT(table.nnz == exp_nnz);
    ASSERT(table.sample_ids == exp_sids);
    ASSERT(table.obs_ids == exp_oids);
    ASSERT(table.sample_indptr == exp_s_indptr);
    ASSERT(table.obs_indptr == exp_o_indptr);

    SUITE_END();
}

void test_biom_get_obs_data() {
    SUITE_START("biom get obs data");

    su::biom table = su::biom("test.biom");
    double exp0[] = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    std::vector<double> exp0_vec = _double_array_to_vector(exp0, 6);
    double exp1[] = {5.0, 1.0, 0.0, 2.0, 3.0, 1.0};
    std::vector<double> exp1_vec = _double_array_to_vector(exp1, 6);
    double exp2[] = {0.0, 0.0, 1.0, 4.0, 0.0, 2.0};
    std::vector<double> exp2_vec = _double_array_to_vector(exp2, 6);
    double exp3[] = {2.0, 1.0, 1.0, 0.0, 0.0, 1.0};
    std::vector<double> exp3_vec = _double_array_to_vector(exp3, 6);
    double exp4[] = {0.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    std::vector<double> exp4_vec = _double_array_to_vector(exp4, 6);
    
    double *out = (double*)malloc(sizeof(double) * 6);
    std::vector<double> obs_vec;

    table.get_obs_data(std::string("GG_OTU_1").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp0_vec));

    table.get_obs_data(std::string("GG_OTU_2").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp1_vec));

    table.get_obs_data(std::string("GG_OTU_3").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp2_vec));

    table.get_obs_data(std::string("GG_OTU_4").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp3_vec));

    table.get_obs_data(std::string("GG_OTU_5").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp4_vec));

    free(out);
    SUITE_END();
}

void test_biom_get_sample_data() {
    SUITE_START("biom get sample data");

    su::biom table = su::biom("test.biom");
    double exp0[] = {0.0, 5.0, 0.0, 2.0, 0.0};
    std::vector<double> exp0_vec = _double_array_to_vector(exp0, 5);
    double exp1[] = {0.0, 1.0, 0.0, 1.0, 1.0};
    std::vector<double> exp1_vec = _double_array_to_vector(exp1, 5);
    double exp2[] = {1.0, 0.0, 1.0, 1.0, 1.0};
    std::vector<double> exp2_vec = _double_array_to_vector(exp2, 5);
    double exp3[] = {0.0, 2.0, 4.0, 0.0, 0.0};
    std::vector<double> exp3_vec = _double_array_to_vector(exp3, 5);
    double exp4[] = {0.0, 3.0, 0.0, 0.0, 0.0};
    std::vector<double> exp4_vec = _double_array_to_vector(exp4, 5);
    double exp5[] = {0.0, 1.0, 2.0, 1.0, 0.0};
    std::vector<double> exp5_vec = _double_array_to_vector(exp5, 5);
    
    double *out = (double*)malloc(sizeof(double) * 5);
    std::vector<double> obs_vec;

    table.get_sample_data(std::string("Sample1").c_str(), out);
    obs_vec = _double_array_to_vector(out, 5);
    ASSERT(vec_almost_equal(obs_vec, exp0_vec));

    table.get_sample_data(std::string("Sample2").c_str(), out);
    obs_vec = _double_array_to_vector(out, 5);
    ASSERT(vec_almost_equal(obs_vec, exp1_vec));

    table.get_sample_data(std::string("Sample3").c_str(), out);
    obs_vec = _double_array_to_vector(out, 5);
    ASSERT(vec_almost_equal(obs_vec, exp2_vec));

    table.get_sample_data(std::string("Sample4").c_str(), out);
    obs_vec = _double_array_to_vector(out, 5);
    ASSERT(vec_almost_equal(obs_vec, exp3_vec));

    table.get_sample_data(std::string("Sample5").c_str(), out);
    obs_vec = _double_array_to_vector(out, 5);
    ASSERT(vec_almost_equal(obs_vec, exp4_vec));

    table.get_sample_data(std::string("Sample6").c_str(), out);
    obs_vec = _double_array_to_vector(out, 5);
    ASSERT(vec_almost_equal(obs_vec, exp5_vec));

    free(out);
    SUITE_END();
}

void test_bptree_leftchild() {
    SUITE_START("test bptree left child");
    su::BPTree tree = su::BPTree("((3,4,(6)5)2,7,((10,100)9)8)1;");

    uint32_t exp[] = {1, 2, 0, 0, 7, 0, 0, 14, 15, 0, 0};
    std::vector<bool> structure = tree.get_structure();

    uint32_t exp_pos = 0;
    for(int i = 0; i < tree.nparens; i++) {
        if(structure[i])
            ASSERT(tree.leftchild(i) == exp[exp_pos++]);
    }
    SUITE_END();
}

void test_bptree_rightchild() {
    SUITE_START("test bptree right child");
    su::BPTree tree = su::BPTree("((3,4,(6)5)2,7,((10,100)9)8)1;");

    uint32_t exp[] = {13, 6, 0, 0, 7, 0, 0, 14, 17, 0, 0};
    std::vector<bool> structure = tree.get_structure();

    uint32_t exp_pos = 0;
    for(int i = 0; i < tree.nparens; i++) {
        if(structure[i])
            ASSERT(tree.rightchild(i) == exp[exp_pos++]);
    }
    SUITE_END();
}

void test_bptree_rightsibling() {
    SUITE_START("test bptree rightsibling");
    su::BPTree tree = su::BPTree("((3,4,(6)5)2,7,((10,100)9)8)1;");

    uint32_t exp[] = {0, 11, 4, 6, 0, 0, 13, 0, 0, 17, 0};
    std::vector<bool> structure = tree.get_structure();
    
    uint32_t exp_pos = 0;
    for(int i = 0; i < tree.nparens; i++) {
        if(structure[i])
            ASSERT(tree.rightsibling(i) == exp[exp_pos++]);
    }
    SUITE_END();
}

void test_propstack_constructor() {
    SUITE_START("test propstack constructor");
    su::PropStack ps = su::PropStack(10);
    // nothing to test directly...
    SUITE_END();
}

void test_propstack_push_and_pop() {
    SUITE_START("test propstack push and pop");
    su::PropStack ps = su::PropStack(10);

    double *vec1 = ps.pop(1);
    double *vec2 = ps.pop(2);
    double *vec3 = ps.pop(3);
    double *vec1_obs;
    double *vec2_obs;
    double *vec3_obs;

    ps.push(1); 
    ps.push(2); 
    ps.push(3);
    
    vec3_obs = ps.pop(4);
    vec2_obs = ps.pop(5);
    vec1_obs = ps.pop(6);

    ASSERT(vec1 == vec1_obs);
    ASSERT(vec2 == vec2_obs);
    ASSERT(vec3 == vec3_obs);
    SUITE_END();
}

void test_propstack_get() {
    SUITE_START("test propstack get");
    su::PropStack ps = su::PropStack(10);

    double *vec1 = ps.pop(1);
    double *vec2 = ps.pop(2);
    double *vec3 = ps.pop(3);
    
    double *vec1_obs = ps.get(1);
    double *vec2_obs = ps.get(2);
    double *vec3_obs = ps.get(3);

    ASSERT(vec1 == vec1_obs);
    ASSERT(vec2 == vec2_obs);
    ASSERT(vec3 == vec3_obs);
    SUITE_END();
}

void test_unifrac_set_proportions() {
    SUITE_START("test unifrac set proportions");
    //                           0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
    //                           ( ( ) ( ( ) ( ) ) ( ( ) ( ) ) )
    su::BPTree tree = su::BPTree("(GG_OTU_1,(GG_OTU_2,GG_OTU_3),(GG_OTU_5,GG_OTU_4));");
    su::biom table = su::biom("test.biom");
    su::PropStack ps = su::PropStack(table.n_samples);
   
    double sample_counts[] = {7, 3, 4, 6, 3, 4}; 
    double *obs = ps.pop(4); // GG_OTU_2
    double exp4[] = {0.714285714286, 0.333333333333, 0.0, 0.333333333333, 1.0, 0.25};
    set_proportions(obs, tree, 4, table, ps, sample_counts);
    for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp4[i]) < 0.000001);

    obs = ps.pop(6); // GG_OTU_3
    double exp6[] = {0.0, 0.0, 0.25, 0.666666666667, 0.0, 0.5};
    set_proportions(obs, tree, 6, table, ps, sample_counts);
    for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp6[i]) < 0.000001);

    obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
    double exp3[] = {0.71428571, 0.33333333, 0.25, 1.0, 1.0, 0.75};
    set_proportions(obs, tree, 3, table, ps, sample_counts);
    for(unsigned int i = 0; i < table.n_samples; i++) 
        ASSERT(fabs(obs[i] - exp3[i]) < 0.000001);
    SUITE_END();
}

void test_unifrac_deconvolute_stripes() {
    SUITE_START("test deconvolute stripes");
    std::vector<double*> stripes;
    double s1[] = {1, 1, 1, 1, 1, 1};
    double s2[] = {2, 2, 2, 2, 2, 2};
    double s3[] = {3, 3, 3, 3, 3, 3};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);

    double exp[6][6] = { {0, 1, 2, 3, 2, 1},
                         {1, 0, 1, 2, 3, 2},
                         {2, 1, 0, 1, 2, 3},
                         {3, 2, 1, 0, 1, 2},
                         {2, 3, 2, 1, 0, 1},
                         {1, 2, 3, 2, 1, 0} };
    double **obs = su::deconvolute_stripes(stripes, 6);
    for(unsigned int i = 0; i < 6; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(exp[i][j] == obs[i][j]);
        }
    }
    free(obs);
    SUITE_END();
}

void test_unnormalized_weighted_unifrac() {
    SUITE_START("test unnormalized weighted unifrac");
    double exp[6][6] = {{ 0.        ,  1.52380952,  2.17857143,  1.9047619 ,  1.14285714, 1.07142857},
                        { 1.52380952,  0.        ,  1.25      ,  2.66666667,  2.66666667, 1.83333333},
                        { 2.17857143,  1.25      ,  0.        ,  2.75      ,  3.25      , 1.75      },
                        { 1.9047619 ,  2.66666667,  2.75      ,  0.        ,  1.33333333, 1.        },
                        { 1.14285714,  2.66666667,  3.25      ,  1.33333333,  0.        , 2.        },
                        { 1.07142857,  1.83333333,  1.75      ,  1.        ,  2.        , 0.        }};
    double **obs;

    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");
    obs = su::unifrac(table, tree, su::weighted_unnormalized);
    for(unsigned int i = 0; i < 6; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(obs[i][j] - exp[i][j]) < 0.000001);
        }
        free(obs[i]);
    }
    SUITE_END();
}

// normalized and unweighted are the same thing short of a conversion to presence/absence
// prior to conversion to proportions.
void test_unweighted_unifrac() {
    SUITE_START("test unweighted unifrac");
    double exp[6][6] = {{ 0.        ,  0.2       ,  0.57142857,  0.6       ,  0.5       , 0.2       },
                        { 0.2       ,  0.        ,  0.42857143,  0.66666667,  0.6       , 0.33333333},
                        { 0.57142857,  0.42857143,  0.        ,  0.71428571,  0.85714286, 0.42857143},
                        { 0.6       ,  0.66666667,  0.71428571,  0.        ,  0.33333333, 0.4       },
                        { 0.5       ,  0.6       ,  0.85714286,  0.33333333,  0.        , 0.6       },
                        { 0.2       ,  0.33333333,  0.42857143,  0.4       ,  0.6       , 0.        }};
    double **obs;
    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");
    obs = su::unifrac(table, tree, su::unweighted);
    for(unsigned int i = 0; i < 6; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(obs[i][j] - exp[i][j]) < 0.000001);
        }
        free(obs[i]);
    }
    SUITE_END();
}

void test_normalized_weighted_unifrac() {
    SUITE_START("test normalized weighted unifrac");
    double exp[6][6] = {{ 0.        ,  0.38095238,  0.58095238,  0.47619048,  0.28571429, 0.26785714},
                        { 0.38095238,  0.        ,  0.33333333,  0.66666667,  0.66666667, 0.45833333},
                        { 0.58095238,  0.33333333,  0.        ,  0.73333333,  0.86666667, 0.46666667},
                        { 0.47619048,  0.66666667,  0.73333333,  0.        ,  0.33333333, 0.25      },
                        { 0.28571429,  0.66666667,  0.86666667,  0.33333333,  0.        , 0.5       },
                        { 0.26785714,  0.45833333,  0.46666667,  0.25      ,  0.5       , 0.        }};
    double **obs;
    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");
    obs = su::unifrac(table, tree, su::weighted_normalized);
    for(unsigned int i = 0; i < 6; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(obs[i][j] - exp[i][j]) < 0.000001);
        }
        free(obs[i]);
    }
    SUITE_END();
}

void test_bptree_shear_simple() {
    SUITE_START("test bptree shear simple");
    su::BPTree tree = su::BPTree("((3:2,4:3,(6:5)5:4)2:1,7:6,((10:9,11:10)9:8)8:7)r");
 
    // simple
    std::unordered_set<std::string> to_keep = {"4", "6", "7", "10", "11"};	
    
    uint32_t exp_nparens = 20;
    std::vector<bool> exp_structure = {true, true, true, false, true, true, false, false, false, true, 
                                       false, true, true, true, false, true, false, false, false, false};
    std::vector<std::string> exp_names = {"r", "2", "4", "", "5", "6", "", "", "", "7", "", "8", "9", "10", "", 
                                          "11", "", "", "", ""};
    std::vector<double> exp_lengths = {0, 1, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0, 10, 0, 0, 0, 0};

    su::BPTree obs = tree.shear(to_keep);
    ASSERT(obs.get_structure() == exp_structure);
    ASSERT(exp_nparens == obs.nparens);
    ASSERT(vec_almost_equal(exp_lengths, obs.lengths));
    ASSERT(obs.names == exp_names);
    SUITE_END();
}

void test_bptree_shear_deep() {
    SUITE_START("test bptree shear deep");
    su::BPTree tree = su::BPTree("((3:2,4:3,(6:5)5:4)2:1,7:6,((10:9,11:10)9:8)8:7)r");
 
    // deep
    std::unordered_set<std::string> to_keep = {"10", "11"};	
    
    uint32_t exp_nparens = 10;
    std::vector<bool> exp_structure = {true, true, true, true, false, true, false, false, false, false};
    std::vector<std::string> exp_names = {"r", "8", "9", "10", "", "11", "", "", "", ""};
    std::vector<double> exp_lengths = {0, 7, 8, 9, 0, 10, 0, 0, 0, 0};

    su::BPTree obs = tree.shear(to_keep);
    ASSERT(exp_nparens == obs.nparens);
    ASSERT(obs.get_structure() == exp_structure);
    ASSERT(vec_almost_equal(exp_lengths, obs.lengths));
    ASSERT(obs.names == exp_names);
    SUITE_END();
}


void test_bptree_collapse_simple() {
    SUITE_START("test bptree collapse simple");
    su::BPTree tree = su::BPTree("((3:2,4:3,(6:5)5:4)2:1,7:6,((10:9,11:10)9:8)8:7)r");
    
    uint32_t exp_nparens = 18;
    std::vector<bool> exp_structure = {true, true, true, false, true, false, true, false, false, 
                                       true, false, true, true, false, true, false, false, false};
    std::vector<std::string> exp_names = {"r", "2", "3", "", "4", "", "6", "", "", "7", "", "9", "10", "", "11", "", "", ""};
    std::vector<double> exp_lengths = {0, 1, 2, 0, 3, 0, 9, 0, 0, 6, 0, 15, 9, 0, 10, 0, 0, 0};

    su::BPTree obs = tree.collapse();

    ASSERT(obs.get_structure() == exp_structure);
    ASSERT(exp_nparens == obs.nparens);
    ASSERT(vec_almost_equal(exp_lengths, obs.lengths));
    ASSERT(obs.names == exp_names);
    SUITE_END();
}

void test_bptree_collapse_edge() {
    SUITE_START("test bptree collapse edge case against root");

    su::BPTree tree = su::BPTree("((a),b)r;");
    su::BPTree exp = su::BPTree("(a,b)r;");
    su::BPTree obs = tree.collapse();
    ASSERT(obs.get_structure() == exp.get_structure());
    ASSERT(obs.names == exp.names);
    ASSERT(vec_almost_equal(obs.lengths, exp.lengths));

	SUITE_END();
}

void test_unifrac_sample_counts() {
    SUITE_START("test unifrac sample counts");
    su::biom table = su::biom("test.biom");
    double* obs = get_sample_counts(table);
    double exp[] = {7, 3, 4, 6, 3, 4};
    for(unsigned int i = 0; i < 6; i++)
        ASSERT(obs[i] == exp[i]);
    SUITE_END();
}

int main(int argc, char** argv) {
    test_bptree_constructor_simple();
    test_bptree_constructor_from_existing();
    test_bptree_constructor_single_descendent();
    test_bptree_constructor_complex();
    test_bptree_constructor_semicolon();
    test_bptree_constructor_edgecases();
    test_bptree_postorder();
    test_bptree_preorder();
    test_bptree_parent();
    test_bptree_leftchild();
    test_bptree_rightchild();
    test_bptree_rightsibling();
    test_bptree_mask();
    test_bptree_shear_simple();
    test_bptree_shear_deep();
    test_bptree_collapse_simple();
    test_bptree_collapse_edge();

    test_biom_constructor();
    test_biom_get_obs_data();
    test_biom_get_sample_data();

    test_propstack_constructor();
    test_propstack_push_and_pop();
    test_propstack_get();

    test_unifrac_set_proportions();
    test_unifrac_deconvolute_stripes();
    test_unweighted_unifrac();
    test_unnormalized_weighted_unifrac();
    test_normalized_weighted_unifrac();
    test_unifrac_sample_counts();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);
  
    printf("\n THE END.\n");
    
    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

#include <iostream>
#include "api.hpp"
#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include "unifrac_internal.hpp"
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

void test_bptree_constructor_quoted_comma() {
    SUITE_START("quoted comma bug");
    su::BPTree tree = su::BPTree("((3,'foo,bar')x,c)r;");
    std::vector<std::string> exp_names = {"r", "x", "3", "", "foo,bar", "", "", "c", "", ""};
    ASSERT(exp_names.size() == tree.names.size());

    for(unsigned int i = 0; i < tree.names.size(); i++) {
        ASSERT(exp_names[i] == tree.names[i]);
    }
    SUITE_END();
}

void test_bptree_constructor_quoted_parens() {
    SUITE_START("quoted parens");
    su::BPTree tree = su::BPTree("((3,'foo(b)ar')x,c)r;");
    std::vector<std::string> exp_names = {"r", "x", "3", "", "foo(b)ar", "", "", "c", "", ""};
    ASSERT(exp_names.size() == tree.names.size());

    for(unsigned int i = 0; i < tree.names.size(); i++) {
        ASSERT(exp_names[i] == tree.names[i]);
    }
    SUITE_END();
}
void test_bptree_postorder() {
    SUITE_START("postorderselect");

    // fig1 from https://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    su::BPTree tree = su::BPTree("((3,4,(6)5)2,7,((10,100)9)8)1;");
    uint32_t exp[] = {2, 4, 7, 6, 1, 11, 15, 17, 14, 13, 0};
    uint32_t obs[tree.nparens / 2];

    for(unsigned int i = 0; i < (tree.nparens / 2); i++)
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

    for(unsigned int i = 0; i < (tree.nparens / 2); i++)
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

    for(int i = 0; i < (int(tree.nparens) - 2); i++)
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

void test_bptree_leftchild() {
    SUITE_START("test bptree left child");
    su::BPTree tree = su::BPTree("((3,4,(6)5)2,7,((10,100)9)8)1;");

    uint32_t exp[] = {1, 2, 0, 0, 7, 0, 0, 14, 15, 0, 0};
    std::vector<bool> structure = tree.get_structure();

    uint32_t exp_pos = 0;
    for(unsigned int i = 0; i < tree.nparens; i++) {
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
    for(unsigned int i = 0; i < tree.nparens; i++) {
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
    for(unsigned int i = 0; i < tree.nparens; i++) {
        if(structure[i])
            ASSERT(tree.rightsibling(i) == exp[exp_pos++]);
    }
    SUITE_END();
}

void test_propstack_constructor() {
    SUITE_START("test propstack constructor");
    su::PropStack<double> ps(10);
    // nothing to test directly...
    SUITE_END();
}

void test_propstack_push_and_pop() {
    SUITE_START("test propstack push and pop");
    su::PropStack<double> ps(10);

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
    su::PropStack<double> ps(10);

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
    su::PropStack<double> ps(table.n_samples);

    double *obs = ps.pop(4); // GG_OTU_2
    double exp4[] = {0.714285714286, 0.333333333333, 0.0, 0.333333333333, 1.0, 0.25};
    set_proportions(obs, tree, 4, table, ps);
    for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp4[i]) < 0.000001);

    obs = ps.pop(6); // GG_OTU_3
    double exp6[] = {0.0, 0.0, 0.25, 0.666666666667, 0.0, 0.5};
    set_proportions(obs, tree, 6, table, ps);
    for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp6[i]) < 0.000001);

    obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
    double exp3[] = {0.71428571, 0.33333333, 0.25, 1.0, 1.0, 0.75};
    set_proportions(obs, tree, 3, table, ps);
    for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp3[i]) < 0.000001);
    SUITE_END();
}

void test_unifrac_set_proportions_range() {
    SUITE_START("test unifrac set proportions range");
    //                           0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
    //                           ( ( ) ( ( ) ( ) ) ( ( ) ( ) ) )
    su::BPTree tree = su::BPTree("(GG_OTU_1,(GG_OTU_2,GG_OTU_3),(GG_OTU_5,GG_OTU_4));");
    su::biom table = su::biom("test.biom");

    const double exp4[] = {0.714285714286, 0.333333333333, 0.0, 0.333333333333, 1.0, 0.25};
    const double exp6[] = {0.0, 0.0, 0.25, 0.666666666667, 0.0, 0.5};
    const double exp3[] = {0.71428571, 0.33333333, 0.25, 1.0, 1.0, 0.75};


    // first the whole table
    {
      su::PropStack<double> ps(table.n_samples);

      double *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, 0, table.n_samples, ps);
      for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, 0, table.n_samples, ps);
      for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, 0, table.n_samples, ps);
      for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp3[i]) < 0.000001);
    }

    // beginning
    {
      su::PropStack<double> ps(3);

      double *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, 0, 3, ps);
      for(unsigned int i = 0; i < 3; i++)
        ASSERT(fabs(obs[i] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, 0, 3, ps);
      for(unsigned int i = 0; i < 3; i++)
        ASSERT(fabs(obs[i] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, 0, 3, ps);
      for(unsigned int i = 0; i < 3; i++)
        ASSERT(fabs(obs[i] - exp3[i]) < 0.000001);
    }

    
    // end
    {
      su::PropStack<double> ps(4);

      double *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, 2, table.n_samples, ps);
      for(unsigned int i = 2; i < table.n_samples; i++)
        ASSERT(fabs(obs[i-2] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, 2, table.n_samples, ps);
      for(unsigned int i = 2; i < table.n_samples; i++)
        ASSERT(fabs(obs[i-2] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, 2, table.n_samples, ps);
      for(unsigned int i = 2; i < table.n_samples; i++)
        ASSERT(fabs(obs[i-2] - exp3[i]) < 0.000001);
    }

    
    // middle
    {
      const unsigned int start = 1;
      const unsigned int end = 4;
      su::PropStack<double> ps(end-start);

      double *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, start, end, ps);
      for(unsigned int i =start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, start, end, ps);
      for(unsigned int i = start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, start, end, ps);
      for(unsigned int i = start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp3[i]) < 0.000001);
    }

    SUITE_END();
}

void test_unifrac_set_proportions_range_float() {
    SUITE_START("test unifrac set proportions range float");
    //                           0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
    //                           ( ( ) ( ( ) ( ) ) ( ( ) ( ) ) )
    su::BPTree tree = su::BPTree("(GG_OTU_1,(GG_OTU_2,GG_OTU_3),(GG_OTU_5,GG_OTU_4));");
    su::biom table = su::biom("test.biom");

    const float exp4[] = {0.714285714286, 0.333333333333, 0.0, 0.333333333333, 1.0, 0.25};
    const float exp6[] = {0.0, 0.0, 0.25, 0.666666666667, 0.0, 0.5};
    const float exp3[] = {0.71428571, 0.33333333, 0.25, 1.0, 1.0, 0.75};

    // just midle
    {
      const unsigned int start = 1;
      const unsigned int end = 4;
      su::PropStack<float> ps(end-start);

      float *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, start, end, ps);
      for(unsigned int i =start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, start, end, ps);
      for(unsigned int i = start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, start, end, ps);
      for(unsigned int i = start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp3[i]) < 0.000001);
    }

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

void test_unifrac_stripes_to_condensed_form_even() {
    SUITE_START("test stripes_to_condensed_form even samples");
    std::vector<double*> stripes;
    double s1[] = {0,  9, 17, 24, 30, 35, 39, 42, 44,  8};
    double s2[] = {1, 10, 18, 25, 31, 36, 40, 43,  7, 16};
    double s3[] = {2, 11, 19, 26, 32, 37, 41,  6, 15, 23};
    double s4[] = {3, 12, 20, 27, 33, 38,  5, 14, 22, 29};
    double s5[] = {4, 13, 21, 28, 34,  4, 13, 21, 28, 34};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[45] = {/* 0, */  0,  1,  2,  3,  4,  5,  6,  7,  8,
                      /* *,  0, */  9, 10, 11, 12, 13, 14, 15, 16,
                      /* *,  *,  0, */ 17, 18, 19, 20, 21, 22, 23,
                      /* *,  *,  *,  0, */ 24, 25, 26, 27, 28, 29,
                      /* *,  *,  *,  *,  0, */ 30, 31, 32, 33, 34,
                      /* *,  *,  *,  *,  *,  0, */ 35, 36, 37, 38,
                      /* *,  *,  *,  *,  *,  *,  0, */ 39, 40, 41,
                      /* *,  *,  *,  *,  *,  *,  *,  0, */ 42, 43,
                      /* *,  *,  *,  *,  *,  *,  *,  *,  0, */ 44};
                      /* *,  *,  *,  *,  *,  *,  *,  *,  *, *,  0 */

    double *obs = (double*)malloc(sizeof(double) * 45);
    su::stripes_to_condensed_form(stripes, 10, obs, 0, 5);
    for(unsigned int i = 0; i < 45; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
    SUITE_END();
}

void test_unifrac_stripes_to_condensed_form_odd() {
    SUITE_START("test stripes_to_condensed_form odd samples");
    std::vector<double*> stripes;
    double s1[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 0};
    double s2[] = {20, 19, 18, 17, 16, 15, 14 ,13, 12, 11, 1};
    double s3[] = {21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 2};
    double s4[] = {40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 3};
    double s5[] = {41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 4};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[55] = {/* 0, */ 1, 20, 21, 40, 41, 47, 33, 29, 11,  0,
                      /* 1,  0, */ 2, 19, 22, 39, 42, 48, 32, 30,  1,
                      /*20,  2,  0, */ 3, 18, 23, 38, 43, 49, 31,  2,
                      /*21, 19,  3,  0, */ 4, 17, 24, 37, 44, 50,  3,
                      /*40, 22, 18,  4,  0, */ 5, 16, 25 ,36, 45 , 4,
                      /*41, 39, 23, 17,  5,  0, */ 6, 15, 26, 35, 46,
                      /*47, 42, 38, 24, 16,  6,  0, */ 7, 14, 27, 34,
                      /*33, 48, 43, 37, 25, 15,  7,  0,*/  8, 13, 28,
                      /*29, 32, 49, 44, 36, 26, 14,  8,  0, */ 9, 12,
                      /*11, 30, 31, 50, 45, 35, 27, 13,  9,  0,*/ 10};
                      /* 0,  1,  2,  3,  4, 46, 34, 28, 12, 10,  0}; */
    double *obs = (double*)malloc(sizeof(double) * 55);
    su::stripes_to_condensed_form(stripes, 11, obs, 0, 5);
    for(unsigned int i = 0; i < 55; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
    SUITE_END();
}

void test_unifrac_stripes_to_condensed_form_odd2() {
    SUITE_START("test stripes_to_condensed_form odd(2) samples");
    std::vector<double*> stripes;
    double s1[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9};
    double s2[] = {18, 17, 16, 15, 14, 13, 12 ,11, 10};
    double s3[] = {19, 20, 21, 22, 23, 24, 25, 26, 27};
    double s4[] = {36, 35, 34, 33, 32, 31, 30, 29, 28};
    double s5[] = {31, 30, 29, 28, 36, 35, 34, 33, 32};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[36] = {/* 0, */ 1, 18, 19, 36, 31, 25, 11,  9,
                      /* 1,  0, */ 2, 17, 20, 35, 30, 26, 10,
                      /*20,  2,  0, */ 3, 16, 21, 34, 29, 27,
                      /*21, 19,  3,  0, */ 4, 15, 22, 33, 28,
                      /*40, 22, 18,  4,  0, */ 5, 14, 23 ,32,
                      /*41, 39, 23, 17,  5,  0, */ 6, 13, 24,
                      /*47, 42, 38, 24, 16,  6,  0, */ 7, 12,
                      /*47, 42, 38, 24, 16,  6,  7, 0, */  8};
                      /* 0,  1,  2,  3,  4, 46, 34, 8,  8, 0}; */
    double *obs = (double*)malloc(sizeof(double) * 36);
    su::stripes_to_condensed_form(stripes, 9, obs, 0, 5);
    for(unsigned int i = 0; i < 36; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
    SUITE_END();
}

class ValidatedMemoryStripes : public su::MemoryStripes {
        private:
           const uint32_t n_stripes;
           mutable std::vector<uint8_t> stripe_status; // 0 new, 1 allocated, 2 deallocated, 3 reallocated, 6 deallocate after rellocation
        public:
           ValidatedMemoryStripes(uint32_t _n_stripes, std::vector<double*> &_stripes) 
           : su::MemoryStripes(_stripes) 
           , n_stripes(_n_stripes)
           , stripe_status(n_stripes)
           {
             for (uint32_t i=0; i<n_stripes; i++) stripe_status[i] = 0;
           }

           virtual const double *get_stripe(uint32_t stripe) const {
              stripe_status[stripe]|=1; 
              return su::MemoryStripes::get_stripe(stripe);
           }
           virtual void release_stripe(uint32_t stripe) const { 
              if (stripe_status[stripe]<2) {
                 stripe_status[stripe]=2;
              } else {
                 stripe_status[stripe]=6;
              }
           }

           bool allInitialized() const {
             bool out = true;
             for (uint32_t i=0; i<n_stripes; i++) out &= (stripe_status[i] != 0);
             return out;
           }

           bool allDealocated() const {
             bool out = true;
             for (uint32_t i=0; i<n_stripes; i++) out &= ((stripe_status[i]&2) == 2);
             return out;
           }

           bool anyRealocated() const {
             bool out = false;
             for (uint32_t i=0; i<n_stripes; i++) out |= (stripe_status[i] >2);
             return out;
           }


};


void test_unifrac_stripes_to_matrix_even() {
    SUITE_START("test stripes_to_matrix even samples");
    std::vector<double*> stripes;
    double s1[] = {0,  9, 17, 24, 30, 35, 39, 42, 44,  8};
    double s2[] = {1, 10, 18, 25, 31, 36, 40, 43,  7, 16};
    double s3[] = {2, 11, 19, 26, 32, 37, 41,  6, 15, 23};
    double s4[] = {3, 12, 20, 27, 33, 38,  5, 14, 22, 29};
    double s5[] = {4, 13, 21, 28, 34,  4, 13, 21, 28, 34};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    // test also double to float conversion
    float exp[100] = {0,  0,  1,  2,  3,  4,  5,  6,  7,  8, 
                      0,  0,  9, 10, 11, 12, 13, 14, 15, 16,  
                      1,  9,  0, 17, 18, 19, 20, 21, 22, 23,
                      2, 10, 17,  0, 24, 25, 26, 27, 28, 29, 
                      3, 11, 18, 24,  0, 30, 31, 32, 33, 34,
                      4, 12, 19, 25, 30,  0, 35, 36, 37, 38,
                      5, 13, 20, 26, 31, 35,  0, 39, 40, 41,
                      6, 14, 21, 27, 32, 36, 39,  0, 42, 43,
                      7, 15, 22, 28, 33, 37, 40, 42,  0, 44,
                      8, 16, 23, 29, 34, 38, 41, 43, 44,  0};
    {
      float *obs = (float*)malloc(sizeof(float) * 100);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix_fp32(vs, 10, 5, obs);
      for(unsigned int i = 0; i < 100; i++) {
        ASSERT(exp[i] == obs[i]);
      }

      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true); 
      ASSERT(vs.anyRealocated() == false);

      free(obs);
    }

    { // small tiles
      float *obs = (float*)malloc(sizeof(float) * 100);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix_fp32(vs, 10, 5, obs, 4);
      for(unsigned int i = 0; i < 100; i++) {
        ASSERT(exp[i] == obs[i]);
      }

      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
    
      free(obs);
    }

    { // large tiles
      float *obs = (float*)malloc(sizeof(float) * 100);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix_fp32(vs, 10, 5, obs, 128);
      for(unsigned int i = 0; i < 100; i++) {
        ASSERT(exp[i] == obs[i]);
      }

      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
    
      free(obs);
    }


    // test also intermediate, 2-step procedure
    double *obsC = (double*)malloc(sizeof(double) * 45);
    su::stripes_to_condensed_form(stripes, 10, obsC, 0, 5);

    float *obs2 = (float*)malloc(sizeof(float) * 100);
    su::condensed_form_to_matrix_fp32(obsC, 10, obs2);

    for(unsigned int i = 0; i < 100; i++) {
        ASSERT(exp[i] == obs2[i]);
    }

    free(obs2);
    free(obsC);
    SUITE_END();
}

void test_unifrac_stripes_to_matrix_odd() {
    SUITE_START("test stripes_to_matrix odd samples");
    std::vector<double*> stripes;
    double s1[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 0};
    double s2[] = {20, 19, 18, 17, 16, 15, 14 ,13, 12, 11, 1};
    double s3[] = {21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 2};
    double s4[] = {40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 3};
    double s5[] = {41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 4};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[121] = { 0,  1, 20, 21, 40, 41, 47, 33, 29, 11,  0, 
                        1,  0,  2, 19, 22, 39, 42, 48, 32, 30,  1,
                       20,  2,  0,  3, 18, 23, 38, 43, 49, 31,  2,
                       21, 19,  3,  0,  4, 17, 24, 37, 44, 50,  3,
                       40, 22, 18,  4,  0,  5, 16, 25 ,36, 45 , 4, 
                       41, 39, 23, 17,  5,  0,  6, 15, 26, 35, 46,
                       47, 42, 38, 24, 16,  6,  0,  7, 14, 27, 34,
                       33, 48, 43, 37, 25, 15,  7,  0,  8, 13, 28, 
                       29, 32, 49, 44, 36, 26, 14,  8,  0,  9, 12,
                       11, 30, 31, 50, 45, 35, 27, 13,  9,  0, 10,
                        0,  1,  2,  3,  4, 46, 34, 28, 12, 10,  0};

    {
      double *obs = (double*)malloc(sizeof(double) * 121);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 11, 5, obs);
      for(unsigned int i = 0; i < 121; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    { // small tiling
      double *obs = (double*)malloc(sizeof(double) * 121);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 11, 5, obs,4);
      for(unsigned int i = 0; i < 121; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    { // large tiling
      double *obs = (double*)malloc(sizeof(double) * 121);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 11, 5, obs,128);
      for(unsigned int i = 0; i < 121; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }


    // test also intermediate, 2-step procedure
    double *obsC = (double*)malloc(sizeof(double) * 55);
    su::stripes_to_condensed_form(stripes, 11, obsC, 0, 5);

    double *obs2 = (double*)malloc(sizeof(double) * 121);
    su::condensed_form_to_matrix(obsC, 11, obs2);
    
    for(unsigned int i = 0; i < 121; i++) {
        ASSERT(exp[i] == obs2[i]);
    }

    free(obs2);
    free(obsC);
    SUITE_END();
}

void test_unifrac_stripes_to_matrix_odd2() {
    SUITE_START("test stripes_to_matrix odd(2) samples");
    std::vector<double*> stripes;
    double s1[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9};
    double s2[] = {18, 17, 16, 15, 14, 13, 12 ,11, 10};
    double s3[] = {19, 20, 21, 22, 23, 24, 25, 26, 27};
    double s4[] = {36, 35, 34, 33, 32, 31, 30, 29, 28};
    double s5[] = {31, 30, 29, 28, 36, 35, 34, 33, 32};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[81] = { 0,  1, 18, 19, 36, 31, 25, 11,  9,
                       1,  0,  2, 17, 20, 35, 30, 26, 10,
                      18,  2,  0,  3, 16, 21, 34, 29, 27,
                      19, 17,  3,  0,  4, 15, 22, 33, 28,
                      36, 20, 16,  4,  0,  5, 14, 23 ,32,
                      31, 35, 21, 15,  5,  0,  6, 13, 24,
                      25, 30, 34, 22, 14,  6,  0,  7, 12,
                      11, 26, 29, 33, 23, 13,  7,  0,  8,
                       9, 10, 27, 28, 32, 24, 12,  8,  0};

    {
      double *obs = (double*)malloc(sizeof(double) * 81);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 9, 5, obs);
      for(unsigned int i = 0; i < 81; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    { // small tile
      double *obs = (double*)malloc(sizeof(double) * 81);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 9, 5, obs,4);
      for(unsigned int i = 0; i < 81; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    { // large tile
      double *obs = (double*)malloc(sizeof(double) * 81);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 9, 5, obs,128);
      for(unsigned int i = 0; i < 81; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    // test also intermediate, 2-step procedure
    double *obsC = (double*)malloc(sizeof(double) * 36);
    su::stripes_to_condensed_form(stripes, 9, obsC, 0, 5);

    double *obs2 = (double*)malloc(sizeof(double) * 81);
    su::condensed_form_to_matrix(obsC, 9, obs2);
    
    for(unsigned int i = 0; i < 81; i++) {
        ASSERT(exp[i] == obs2[i]);
    }

    free(obs2);
    free(obsC);
    SUITE_END();
}


void test_unnormalized_weighted_unifrac() {
    SUITE_START("test unnormalized weighted unifrac");

    std::vector<std::thread> threads(1);
    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");

    std::vector<double*> exp;
    double stride1[] = {1.52380952, 1.25, 2.75, 1.33333333, 2., 1.07142857};
    double stride2[] = {2.17857143, 2.66666667, 3.25, 1.0, 1.14285714, 1.83333333};
    double stride3[] = {1.9047619, 2.66666667, 1.75, 1.9047619, 2.66666667, 1.75};
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = false;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::weighted_unnormalized,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(threads),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
    }
    SUITE_END();
}

void test_generalized_unifrac() {
    SUITE_START("test generalized unifrac");

    std::vector<std::thread> threads(1);
    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");

    // weighted normalized unifrac as computed above
    std::vector<double*> w_exp;
    double w_stride1[] = {0.38095238, 0.33333333, 0.73333333, 0.33333333, 0.5, 0.26785714};
    double w_stride2[] = {0.58095238, 0.66666667, 0.86666667, 0.25, 0.28571429, 0.45833333};
    double w_stride3[] = {0.47619048, 0.66666667, 0.46666667, 0.47619048, 0.66666667, 0.46666667};
    w_exp.push_back(w_stride1);
    w_exp.push_back(w_stride2);
    w_exp.push_back(w_stride3);
    std::vector<double*> w_strides = su::make_strides(6);
    std::vector<double*> w_strides_total = su::make_strides(6);
    su::task_parameters w_task_p;
    w_task_p.start = 0; w_task_p.stop = 3; w_task_p.tid = 0; w_task_p.n_samples = 6; w_task_p.bypass_tips = false;
    w_task_p.g_unifrac_alpha = 1.0;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(w_task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::generalized,
                        false,
                        std::ref(w_strides),
                        std::ref(w_strides_total),
                        std::ref(threads),
                        std::ref(tasks));

    // as computed by GUniFrac v1.0
    //          Sample1   Sample2   Sample3   Sample4   Sample5   Sample6
    //Sample1 0.0000000 0.4408392 0.6886965 0.7060606 0.5833333 0.3278410
    //Sample2 0.4408392 0.0000000 0.5102041 0.7500000 0.8000000 0.5208125
    //Sample3 0.6886965 0.5102041 0.0000000 0.8649351 0.9428571 0.5952381
    //Sample4 0.7060606 0.7500000 0.8649351 0.0000000 0.5000000 0.4857143
    //Sample5 0.5833333 0.8000000 0.9428571 0.5000000 0.0000000 0.7485714
    //Sample6 0.3278410 0.5208125 0.5952381 0.4857143 0.7485714 0.0000000
    std::vector<double*> d0_exp;
    double d0_stride1[] = {0.4408392, 0.5102041, 0.8649351, 0.5000000, 0.7485714, 0.3278410};
    double d0_stride2[] = {0.6886965, 0.7500000, 0.9428571, 0.4857143, 0.5833333, 0.5208125};
    double d0_stride3[] = {0.7060606, 0.8000000, 0.5952381, 0.7060606, 0.8000000, 0.5952381};
    d0_exp.push_back(d0_stride1);
    d0_exp.push_back(d0_stride2);
    d0_exp.push_back(d0_stride3);
    std::vector<double*> d0_strides = su::make_strides(6);
    std::vector<double*> d0_strides_total = su::make_strides(6);
    su::task_parameters d0_task_p;
    d0_task_p.start = 0; d0_task_p.stop = 3; d0_task_p.tid = 0; d0_task_p.n_samples = 6; d0_task_p.bypass_tips = false;
    d0_task_p.g_unifrac_alpha = 0.0;

    tasks.clear();
    tasks.push_back(d0_task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::generalized,
                        false,
                        std::ref(d0_strides),
                        std::ref(d0_strides_total),
                        std::ref(threads),
                        std::ref(tasks));

    // as computed by GUniFrac v1.0
    //          Sample1   Sample2   Sample3   Sample4   Sample5   Sample6
    //Sample1 0.0000000 0.4040518 0.6285560 0.5869439 0.4082483 0.2995673
    //Sample2 0.4040518 0.0000000 0.4160597 0.7071068 0.7302479 0.4860856
    //Sample3 0.6285560 0.4160597 0.0000000 0.8005220 0.9073159 0.5218198
    //Sample4 0.5869439 0.7071068 0.8005220 0.0000000 0.4117216 0.3485667
    //Sample5 0.4082483 0.7302479 0.9073159 0.4117216 0.0000000 0.6188282
    //Sample6 0.2995673 0.4860856 0.5218198 0.3485667 0.6188282 0.0000000
    std::vector<double*> d05_exp;
    double d05_stride1[] = {0.4040518, 0.4160597, 0.8005220, 0.4117216, 0.6188282, 0.2995673};
    double d05_stride2[] = {0.6285560, 0.7071068, 0.9073159, 0.3485667, 0.4082483, 0.4860856};
    double d05_stride3[] = {0.5869439, 0.7302479, 0.5218198, 0.5869439, 0.7302479, 0.5218198};
    d05_exp.push_back(d05_stride1);
    d05_exp.push_back(d05_stride2);
    d05_exp.push_back(d05_stride3);
    std::vector<double*> d05_strides = su::make_strides(6);
    std::vector<double*> d05_strides_total = su::make_strides(6);
    su::task_parameters d05_task_p;
    d05_task_p.start = 0; d05_task_p.stop = 3; d05_task_p.tid = 0; d05_task_p.n_samples = 6; d05_task_p.bypass_tips = false;
    d05_task_p.g_unifrac_alpha = 0.5;

    tasks.clear();
    tasks.push_back(d05_task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::generalized,
                        false,
                        std::ref(d05_strides),
                        std::ref(d05_strides_total),
                        std::ref(threads),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(w_strides[i][j] - w_exp[i][j]) < 0.000001);
            ASSERT(fabs(d0_strides[i][j] - d0_exp[i][j]) < 0.000001);
            ASSERT(fabs(d05_strides[i][j] - d05_exp[i][j]) < 0.000001);
        }
        free(w_strides[i]);
        free(d0_strides[i]);
        free(d05_strides[i]);
    }
    SUITE_END();
}

void test_vaw_unifrac_weighted_normalized() {
    SUITE_START("test vaw weighted normalized unifrac");

    std::vector<std::thread> threads(1);
    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");

    // as computed by GUniFrac, the original implementation of VAW-UniFrac
    // could not be found.
    //          Sample1   Sample2   Sample3   Sample4   Sample5   Sample6
    //Sample1 0.0000000 0.4086040 0.6240185 0.4639481 0.2857143 0.2766318
    //Sample2 0.4086040 0.0000000 0.3798594 0.6884992 0.6807616 0.4735781
    //Sample3 0.6240185 0.3798594 0.0000000 0.7713254 0.8812897 0.5047114
    //Sample4 0.4639481 0.6884992 0.7713254 0.0000000 0.6666667 0.2709298
    //Sample5 0.2857143 0.6807616 0.8812897 0.6666667 0.0000000 0.4735991
    //Sample6 0.2766318 0.4735781 0.5047114 0.2709298 0.4735991 0.0000000
    // weighted normalized unifrac as computed above

    std::vector<double*> w_exp;
    double w_stride1[] = {0.4086040, 0.3798594, 0.7713254, 0.6666667, 0.4735991, 0.2766318};
    double w_stride2[] = {0.6240185, 0.6884992, 0.8812897, 0.2709298, 0.2857143, 0.4735781};
    double w_stride3[] = {0.4639481, 0.6807616, 0.5047114, 0.4639481, 0.6807616, 0.5047114};
    w_exp.push_back(w_stride1);
    w_exp.push_back(w_stride2);
    w_exp.push_back(w_stride3);
    std::vector<double*> w_strides = su::make_strides(6);
    std::vector<double*> w_strides_total = su::make_strides(6);
    su::task_parameters w_task_p;
    w_task_p.start = 0; w_task_p.stop = 3; w_task_p.tid = 0; w_task_p.n_samples = 6; w_task_p.bypass_tips = false;
    w_task_p.g_unifrac_alpha = 1.0;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(w_task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::weighted_normalized,
                        true,
                        std::ref(w_strides),
                        std::ref(w_strides_total),
                        std::ref(threads),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(w_strides[i][j] - w_exp[i][j]) < 0.000001);
        }
        free(w_strides[i]);
    }
    SUITE_END();
}


void test_make_strides() {
    SUITE_START("test make stripes");
    std::vector<double*> exp;
    double stride[] = {0., 0., 0.};
    exp.push_back(stride);
    exp.push_back(stride);
    exp.push_back(stride);

    std::vector<double*> obs = su::make_strides(3);
    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(obs[i][j] - exp[i][j]) < 0.000001);
        }
        free(obs[i]);
    }
}

void test_faith_pd() {
    SUITE_START("test faith PD");

    // Note this tree is binary (opposed to example below)
    su::BPTree tree = su::BPTree("((GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1):2,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");

    // make vector of expectations from faith PD
    double exp[6] = {6., 7., 8., 5., 4., 7.};

    // run faith PD to get obs
    double obs[6] = {0, 0, 0, 0, 0, 0};

    su::faith_pd(table, tree, obs);

    // ASSERT that results = expectation
    for (unsigned int i = 0; i < 6; i++){
        ASSERT(fabs(exp[i]-obs[i]) < 0.000001)
    }
    SUITE_END();
}

void test_faith_pd_shear(){
    SUITE_START("test faith PD extra OTUs in tree");

    su::BPTree tree = su::BPTree("((GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1,GG_OTU_ex:9):1):2,(GG_OTU_5:1,GG_OTU_4:1,GG_OTU_ex2:12):1);");
    su::biom table = su::biom("test.biom");

    // make vector of expectations from faith PD
    double exp[6] = {6., 7., 8., 5., 4., 7.};

    // run faith PD to get obs
    double obs[6] = {0, 0, 0, 0, 0, 0};

    std::unordered_set<std::string> to_keep(table.obs_ids.begin(),           \
                                            table.obs_ids.end());            \
    su::BPTree tree_sheared = tree.shear(to_keep).collapse();
    su::faith_pd(table, tree_sheared, obs);

    // ASSERT that results = expectation
    for (unsigned int i = 0; i < 6; i++){
        ASSERT(fabs(exp[i]-obs[i]) < 0.000001)
    }
    SUITE_END();
}

void test_unweighted_unifrac() {
    SUITE_START("test unweighted unifrac");
    std::vector<std::thread> threads(1);
    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");

    std::vector<double*> exp;
    double stride1[] = {0.2, 0.42857143, 0.71428571, 0.33333333, 0.6, 0.2};
    double stride2[] = {0.57142857, 0.66666667, 0.85714286, 0.4, 0.5, 0.33333333};
    double stride3[] = {0.6, 0.6, 0.42857143, 0.6, 0.6, 0.42857143};
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = false;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::unweighted,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(threads),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
    }
    SUITE_END();
}

void test_unweighted_unifrac_fast() {
    SUITE_START("test unweighted unifrac no tips");
    std::vector<std::thread> threads(1);
    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");

    std::vector<double*> exp;
    double stride1[] = {0., 0., 0.5, 0., 0.5, 0.};
    double stride2[] = {0., 0.5, 0.5, 0.5, 0.5, 0.};
    double stride3[] = {0.5, 0.5, 0., 0.5, 0.5, 0.};
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = true;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::unweighted,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(threads),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
    }
    SUITE_END();
}

void test_normalized_weighted_unifrac() {
    SUITE_START("test normalized weighted unifrac");
    std::vector<std::thread> threads(1);
    su::BPTree tree = su::BPTree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table = su::biom("test.biom");

    std::vector<double*> exp;
    double stride1[] = {0.38095238, 0.33333333, 0.73333333, 0.33333333, 0.5, 0.26785714};
    double stride2[] = {0.58095238, 0.66666667, 0.86666667, 0.25, 0.28571429, 0.45833333};
    double stride3[] = {0.47619048, 0.66666667, 0.46666667, 0.47619048, 0.66666667, 0.46666667};
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = false;


    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::weighted_normalized,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(threads),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
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

void test_test_table_ids_are_subset_of_tree() {
    SUITE_START("test test_table_ids_are_subset_of_tree");

    su::BPTree tree = su::BPTree("(a:1,b:2)r;");
    su::biom table = su::biom("test.biom");
    std::string expected = "GG_OTU_1";
    std::string observed = su::test_table_ids_are_subset_of_tree(table, tree);
    ASSERT(observed == expected);

    su::BPTree tree2 = su::BPTree("(GG_OTU_1,GG_OTU_5,GG_OTU_6,GG_OTU_2,GG_OTU_3,GG_OTU_4);");
    su::biom table2 = su::biom("test.biom");
    expected = "";
    observed = su::test_table_ids_are_subset_of_tree(table2, tree2);
    ASSERT(observed == expected);
    SUITE_END();
}


void test_bptree_get_tip_names() {
    SUITE_START("test bptree get_tip_names");
    su::BPTree tree = su::BPTree("((a:2,b:3,(c:5)d:4)e:1,f:6,((g:9,h:10)i:8)j:7)r");

    std::unordered_set<std::string> expected = {"a", "b", "c", "f", "g", "h"};
    std::unordered_set<std::string> observed = tree.get_tip_names();
    ASSERT(observed == expected);
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
    double* obs = table.sample_counts;
    double exp[] = {7, 3, 4, 6, 3, 4};
    for(unsigned int i = 0; i < 6; i++)
        ASSERT(obs[i] == exp[i]);
    SUITE_END();
}

void test_set_tasks() {
    SUITE_START("test set tasks");
    std::vector<su::task_parameters> obs(1);
    std::vector<su::task_parameters> exp(1);

    exp[0].g_unifrac_alpha = 1.0;
    exp[0].n_samples = 100;
    exp[0].bypass_tips = false;
    exp[0].start = 0;
    exp[0].stop = 100;
    exp[0].tid = 0;

    set_tasks(obs, 1.0, 100, 0, 100, false, 1);
    ASSERT(obs[0].g_unifrac_alpha == exp[0].g_unifrac_alpha);
    ASSERT(obs[0].n_samples == exp[0].n_samples);
    ASSERT(obs[0].start == exp[0].start);
    ASSERT(obs[0].stop == exp[0].stop);
    ASSERT(obs[0].tid == exp[0].tid);

    std::vector<su::task_parameters> obs2(2);
    std::vector<su::task_parameters> exp2(2);

    exp2[0].g_unifrac_alpha = 1.0;
    exp2[0].n_samples = 100;
    exp2[0].bypass_tips = false;
    exp2[0].start = 0;
    exp2[0].stop = 50;
    exp2[0].tid = 0;
    exp2[1].g_unifrac_alpha = 1.0;
    exp2[1].n_samples = 100;
    exp2[1].bypass_tips = false;
    exp2[1].start = 50;
    exp2[1].stop = 100;
    exp2[1].tid = 1;

    set_tasks(obs2, 1.0, 100, 0, 100, false, 2);
    for(unsigned int i=0; i < 2; i++) {
        ASSERT(obs2[i].g_unifrac_alpha == exp2[i].g_unifrac_alpha);
        ASSERT(obs2[i].n_samples == exp2[i].n_samples);
        ASSERT(obs2[i].start == exp2[i].start);
        ASSERT(obs2[i].stop == exp2[i].stop);
        ASSERT(obs2[i].tid == exp2[i].tid);
    }

    std::vector<su::task_parameters> obs3(3);
    std::vector<su::task_parameters> exp3(3);

    exp3[0].g_unifrac_alpha = 1.0;
    exp3[0].n_samples = 100;
    exp3[0].bypass_tips = false;
    exp3[0].start = 25;
    exp3[0].stop = 50;
    exp3[0].tid = 0;
    exp3[1].g_unifrac_alpha = 1.0;
    exp3[1].n_samples = 100;
    exp3[1].bypass_tips = false;
    exp3[1].start = 50;
    exp3[1].stop = 75;
    exp3[1].tid = 1;
    exp3[2].g_unifrac_alpha = 1.0;
    exp3[2].n_samples = 100;
    exp3[2].bypass_tips = false;
    exp3[2].start = 75;
    exp3[2].stop = 100;
    exp3[2].tid = 2;

    set_tasks(obs3, 1.0, 100, 25, 100, false, 3);
    for(unsigned int i=0; i < 3; i++) {
        ASSERT(obs3[i].g_unifrac_alpha == exp3[i].g_unifrac_alpha);
        ASSERT(obs3[i].n_samples == exp3[i].n_samples);
        ASSERT(obs3[i].start == exp3[i].start);
        ASSERT(obs3[i].stop == exp3[i].stop);
        ASSERT(obs3[i].tid == exp3[i].tid);
    }

    std::vector<su::task_parameters> obs4(3);
    std::vector<su::task_parameters> exp4(3);

    exp4[0].g_unifrac_alpha = 1.0;
    exp4[0].n_samples = 100;
    exp4[0].bypass_tips = false;
    exp4[0].start = 26;
    exp4[0].stop = 51;
    exp4[0].tid = 0;
    exp4[1].g_unifrac_alpha = 1.0;
    exp4[1].n_samples = 100;
    exp4[1].bypass_tips = false;
    exp4[1].start = 51;
    exp4[1].stop = 76;
    exp4[1].tid = 1;
    exp4[2].g_unifrac_alpha = 1.0;
    exp4[2].n_samples = 100;
    exp4[2].bypass_tips = false;
    exp4[2].start = 76;
    exp4[2].stop = 100;
    exp4[2].tid = 2;

    set_tasks(obs4, 1.0, 100, 26, 100, false, 3);
    for(unsigned int i=0; i < 3; i++) {
        ASSERT(obs4[i].g_unifrac_alpha == exp4[i].g_unifrac_alpha);
        ASSERT(obs4[i].n_samples == exp4[i].n_samples);
        ASSERT(obs4[i].start == exp4[i].start);
        ASSERT(obs4[i].stop == exp4[i].stop);
        ASSERT(obs4[i].tid == exp4[i].tid);
    }

    // set_tasks boundary bug
    std::vector<su::task_parameters> obs16(16);
    std::vector<su::task_parameters> exp16(16);
    set_tasks(obs16, 1.0, 9511, 0, 0, false, 16);
    exp16[15].start = 4459;
    exp16[15].stop = 4756;
    ASSERT(obs16[15].start == exp16[15].start);
    ASSERT(obs16[15].stop == exp16[15].stop);
    SUITE_END();
}

void test_bptree_constructor_newline_bug() {
    SUITE_START("test bptree constructor newline bug");
    su::BPTree tree = su::BPTree("((362be41f31fd26be95ae43a8769b91c0:0.116350803,(a16679d5a10caa9753f171977552d920:0.105836235,((a7acc2abb505c3ee177a12e514d3b994:0.008268754,(4e22aa3508b98813f52e1a12ffdb74ad:0.03144211,8139c4ac825dae48454fb4800fb87896:0.043622957)0.923:0.046588301)0.997:0.120902074,((2d3df7387323e2edcbbfcb6e56a02710:0.031543994,3f6752aabcc291b67a063fb6492fd107:0.091571442)0.759:0.016335166,((d599ebe277afb0dfd4ad3c2176afc50e:5e-09,84d0affc7243c7d6261f3a7d680b873f:0.010245188)0.883:0.048993011,51121722488d0c3da1388d1b117cd239:0.119447926)0.763:0.035660204)0.921:0.058191474)0.776:0.02854575)0.657:0.052060833)0.658:0.032547569,(99647b51f775c8ddde8ed36a7d60dbcd:0.173334268,(f18a9c8112372e2916a66a9778f3741b:0.194813398,(5833416522de0cca717a1abf720079ac:5e-09,(2bf1067d2cd4f09671e3ebe5500205ca:0.031692682,(b32621bcd86cb99e846d8f6fee7c9ab8:0.031330707,1016319c25196d73bdb3096d86a9df2f:5e-09)0.058:0.01028612)0.849:0.010284866)0.791:0.041353384)0.922:0.109470534):0.022169824000000005)root;\n\n");
    SUITE_END();
}

int main(int argc, char** argv) {
    test_bptree_constructor_simple();
    test_bptree_constructor_newline_bug();
    test_bptree_constructor_from_existing();
    test_bptree_constructor_single_descendent();
    test_bptree_constructor_complex();
    test_bptree_constructor_semicolon();
    test_bptree_constructor_edgecases();
    test_bptree_constructor_quoted_comma();
    test_bptree_constructor_quoted_parens();
    test_bptree_postorder();
    test_bptree_preorder();
    test_bptree_parent();
    test_bptree_leftchild();
    test_bptree_rightchild();
    test_bptree_rightsibling();
    test_bptree_get_tip_names();
    test_bptree_mask();
    test_bptree_shear_simple();
    test_bptree_shear_deep();
    test_bptree_collapse_simple();
    test_bptree_collapse_edge();

    test_biom_constructor();
    test_biom_get_obs_data();

    test_propstack_constructor();
    test_propstack_push_and_pop();
    test_propstack_get();

    test_unifrac_set_proportions();
    test_unifrac_set_proportions_range();
    test_unifrac_set_proportions_range_float();
    test_unifrac_deconvolute_stripes();
    test_unifrac_stripes_to_condensed_form_even();
    test_unifrac_stripes_to_condensed_form_odd();
    test_unifrac_stripes_to_condensed_form_odd2();
    test_unifrac_stripes_to_matrix_even();
    test_unifrac_stripes_to_matrix_odd();
    test_unifrac_stripes_to_matrix_odd2();
    test_unweighted_unifrac();
    test_unweighted_unifrac_fast();
    test_unnormalized_weighted_unifrac();
    test_normalized_weighted_unifrac();
    test_generalized_unifrac();
    test_vaw_unifrac_weighted_normalized();
    test_unifrac_sample_counts();
    test_set_tasks();
    test_test_table_ids_are_subset_of_tree();

    test_faith_pd();
    test_faith_pd_shear();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);

    printf("\n THE END.\n");

    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

#include <iostream>
#include "api.hpp"
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


//void test_write_mat() {
//    SUITE_START("test write mat_t");
//    SUITE_END();
//}
//
//void test_read_mat() {
//    SUITE_START("test read mat_t");
//    SUITE_END();
//}
//
partial_mat_t* make_test_pm() {
    partial_mat_t* pm = (partial_mat_t*)malloc(sizeof(partial_mat_t));
    pm->n_samples = 6;
    pm->sample_ids = (char**)malloc(sizeof(char*) * 6);
    pm->sample_ids[0] = (char*)malloc(sizeof(char) * 2);
    pm->sample_ids[0][0] = 'A'; pm->sample_ids[0][1] = '\0';
    pm->sample_ids[1] = (char*)malloc(sizeof(char) * 2);
    pm->sample_ids[1][0] = 'B'; pm->sample_ids[1][1] = '\0';
    pm->sample_ids[2] = (char*)malloc(sizeof(char) * 3);
    pm->sample_ids[2][0] = 'C'; pm->sample_ids[2][1] = 'x'; pm->sample_ids[2][2] = '\0';
    pm->sample_ids[3] = (char*)malloc(sizeof(char) * 2);
    pm->sample_ids[3][0] = 'D'; pm->sample_ids[3][1] = '\0';
    pm->sample_ids[4] = (char*)malloc(sizeof(char) * 2);
    pm->sample_ids[4][0] = 'E'; pm->sample_ids[4][1] = '\0';
    pm->sample_ids[5] = (char*)malloc(sizeof(char) * 2);
    pm->sample_ids[5][0] = 'F'; pm->sample_ids[5][1] = '\0';
    pm->stripes = (double**)malloc(sizeof(double*) * 3);
    pm->stripes[0] = (double*)malloc(sizeof(double) * 6);
    pm->stripes[0][0] = 1; pm->stripes[0][1] = 2; pm->stripes[0][2] = 3; pm->stripes[0][3] = 4; pm->stripes[0][4] = 5; pm->stripes[0][5] = 6;
    pm->stripes[1] = (double*)malloc(sizeof(double) * 6);
    pm->stripes[1][0] = 7; pm->stripes[1][1] = 8; pm->stripes[1][2] = 9; pm->stripes[1][3] = 10; pm->stripes[1][4] = 11; pm->stripes[1][5] = 12;
    pm->stripes[2] = (double*)malloc(sizeof(double) * 6);
    pm->stripes[2][0] = 13; pm->stripes[2][1] = 14; pm->stripes[2][2] = 15; pm->stripes[2][3] = 16; pm->stripes[2][4] = 17; pm->stripes[2][5] = 18;

    return pm;
}

mat_t* mat_three_rep() {
    mat_t* res = (mat_t*)malloc(sizeof(mat_t));
    res->n_samples = 6;
    res->cf_size = 15;
    res->is_upper_triangle = true;
    res->condensed_form = (double*)malloc(sizeof(double) * 15);
    // using second half of third stripe. the last stripe when operating on even numbers of samples is normally redundant with the first half,
    // but that was more annoying in to write up in the tests.
    res->condensed_form[0] = 1;  res->condensed_form[1] = 7;  res->condensed_form[2] = 16; res->condensed_form[3] = 11; res->condensed_form[4] = 6; 
                                 res->condensed_form[5] = 2;  res->condensed_form[6] = 8;  res->condensed_form[7] = 17;  res->condensed_form[8] = 12;  
                                                              res->condensed_form[9] = 3;  res->condensed_form[10] = 9; res->condensed_form[11] = 18; 
                                                                                           res->condensed_form[12] = 4; res->condensed_form[13] = 10;  
                                                                                                                        res->condensed_form[14] = 5;
    res->sample_ids = (char**)malloc(sizeof(char*) * 6);
    res->sample_ids[0] = (char*)malloc(sizeof(char) * 2);
    res->sample_ids[0][0] = 'A'; res->sample_ids[0][1] = '\0';
    res->sample_ids[1] = (char*)malloc(sizeof(char) * 2);
    res->sample_ids[1][0] = 'B'; res->sample_ids[1][1] = '\0';
    res->sample_ids[2] = (char*)malloc(sizeof(char) * 3);
    res->sample_ids[2][0] = 'C'; res->sample_ids[2][1] = 'x'; res->sample_ids[2][2] = '\0';
    res->sample_ids[3] = (char*)malloc(sizeof(char) * 2);
    res->sample_ids[3][0] = 'D'; res->sample_ids[3][1] = '\0';
    res->sample_ids[4] = (char*)malloc(sizeof(char) * 2);
    res->sample_ids[4][0] = 'E'; res->sample_ids[4][1] = '\0';
    res->sample_ids[5] = (char*)malloc(sizeof(char) * 2);
    res->sample_ids[5][0] = 'F'; res->sample_ids[5][1] = '\0';

    return res;
}

void test_read_write_partial_mat() {
    SUITE_START("test read/write partial_mat_t");

    partial_mat_t* pm = make_test_pm();
    pm->stripe_start = 0;
    pm->stripe_stop = 3;
    pm->stripe_total = 3;
    pm->is_upper_triangle = true;
    
    io_status err = write_partial("/tmp/ssu_io.dat", pm);
    ASSERT(err == write_okay);

    partial_mat_t *obs = NULL;
    err = read_partial("/tmp/ssu_io.dat", &obs);
    
    ASSERT(err == read_okay);
    ASSERT(obs->n_samples == 6);
    ASSERT(obs->stripe_start == 0);
    ASSERT(obs->stripe_stop == 3);
    ASSERT(obs->stripe_total == 3);
    ASSERT(strcmp(obs->sample_ids[0], "A") == 0);
    ASSERT(strcmp(obs->sample_ids[1], "B") == 0);
    ASSERT(strcmp(obs->sample_ids[2], "Cx") == 0);
    ASSERT(strcmp(obs->sample_ids[3], "D") == 0);
    ASSERT(strcmp(obs->sample_ids[4], "E") == 0);
    ASSERT(strcmp(obs->sample_ids[5], "F") == 0);

    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 6; j++) {
            ASSERT(obs->stripes[i][j] == ((i * 6) + j + 1));
        }
    }
    SUITE_END();
}

void test_merge_partial_mat() {
    SUITE_START("test merge partial_mat_t");

    // the easy test
    partial_mat_t* pm1 = (partial_mat_t*)malloc(sizeof(partial_mat_t));
    pm1->n_samples = 6;
    pm1->sample_ids = (char**)malloc(sizeof(char*) * 6);
    pm1->sample_ids[0] = (char*)malloc(sizeof(char) * 2);
    pm1->sample_ids[0][0] = 'A'; pm1->sample_ids[0][1] = '\0';
    pm1->sample_ids[1] = (char*)malloc(sizeof(char) * 2);
    pm1->sample_ids[1][0] = 'B'; pm1->sample_ids[1][1] = '\0';
    pm1->sample_ids[2] = (char*)malloc(sizeof(char) * 3);
    pm1->sample_ids[2][0] = 'C'; pm1->sample_ids[2][1] = 'x'; pm1->sample_ids[2][2] = '\0';
    pm1->sample_ids[3] = (char*)malloc(sizeof(char) * 2);
    pm1->sample_ids[3][0] = 'D'; pm1->sample_ids[3][1] = '\0';
    pm1->sample_ids[4] = (char*)malloc(sizeof(char) * 2);
    pm1->sample_ids[4][0] = 'E'; pm1->sample_ids[4][1] = '\0';
    pm1->sample_ids[5] = (char*)malloc(sizeof(char) * 2);
    pm1->sample_ids[5][0] = 'F'; pm1->sample_ids[5][1] = '\0';
    pm1->stripes = (double**)malloc(sizeof(double*) * 2);
    pm1->stripes[0] = (double*)malloc(sizeof(double) * 6);
    pm1->stripes[0][0] = 1; pm1->stripes[0][1] = 2; pm1->stripes[0][2] = 3; pm1->stripes[0][3] = 4; pm1->stripes[0][4] = 5; pm1->stripes[0][5] = 6;
    pm1->stripes[1] = (double*)malloc(sizeof(double) * 6);
    pm1->stripes[1][0] = 7; pm1->stripes[1][1] = 8; pm1->stripes[1][2] = 9; pm1->stripes[1][3] = 10; pm1->stripes[1][4] = 11; pm1->stripes[1][5] = 12;
    pm1->stripe_start = 0;
    pm1->stripe_stop = 2;
    pm1->stripe_total = 3;
    pm1->is_upper_triangle = true;
    
    partial_mat_t* pm2 = (partial_mat_t*)malloc(sizeof(partial_mat_t));
    pm2->n_samples = 6;
    pm2->sample_ids = (char**)malloc(sizeof(char*) * 6);
    pm2->sample_ids[0] = (char*)malloc(sizeof(char) * 2);
    pm2->sample_ids[0][0] = 'A'; pm2->sample_ids[0][1] = '\0';
    pm2->sample_ids[1] = (char*)malloc(sizeof(char) * 2);
    pm2->sample_ids[1][0] = 'B'; pm2->sample_ids[1][1] = '\0';
    pm2->sample_ids[2] = (char*)malloc(sizeof(char) * 3);
    pm2->sample_ids[2][0] = 'C'; pm2->sample_ids[2][1] = 'x'; pm2->sample_ids[2][2] = '\0';
    pm2->sample_ids[3] = (char*)malloc(sizeof(char) * 2);
    pm2->sample_ids[3][0] = 'D'; pm2->sample_ids[3][1] = '\0';
    pm2->sample_ids[4] = (char*)malloc(sizeof(char) * 2);
    pm2->sample_ids[4][0] = 'E'; pm2->sample_ids[4][1] = '\0';
    pm2->sample_ids[5] = (char*)malloc(sizeof(char) * 2);
    pm2->sample_ids[5][0] = 'F'; pm2->sample_ids[5][1] = '\0';
    pm2->stripes = (double**)malloc(sizeof(double*) * 1);
    pm2->stripes[0] = (double*)malloc(sizeof(double) * 6);
    pm2->stripes[0][0] = 13; pm2->stripes[0][1] = 14; pm2->stripes[0][2] = 15; pm2->stripes[0][3] = 16; pm2->stripes[0][4] = 17; pm2->stripes[0][5] = 18;
    pm2->stripe_start = 2;
    pm2->stripe_stop = 3;
    pm2->stripe_total = 3;
    pm2->is_upper_triangle = true;

    mat_t* exp = mat_three_rep();

    partial_mat_t* pms[2];
    pms[0] = pm1;
    pms[1] = pm2;

    mat_t* obs = NULL;
    merge_status err = merge_partial(pms, 2, 1, &obs);
    ASSERT(err == merge_okay);
    ASSERT(obs->cf_size == exp->cf_size);
    ASSERT(obs->n_samples == exp->n_samples);
    ASSERT(obs->is_upper_triangle == exp->is_upper_triangle);
    for(int i = 0; i < obs->cf_size; i++) {
        ASSERT(obs->condensed_form[i] == exp->condensed_form[i]);
    }
    for(int i = 0; i < obs->n_samples; i++)
        ASSERT(strcmp(obs->sample_ids[i], exp->sample_ids[i]) == 0);
  
    // out of order test

    pms[0] = pm2;
    pms[1] = pm1;

    obs = NULL;
    err = merge_partial(pms, 2, 1, &obs);
    ASSERT(err == merge_okay);
    ASSERT(obs->cf_size == exp->cf_size);
    ASSERT(obs->n_samples == exp->n_samples);
    ASSERT(obs->is_upper_triangle == exp->is_upper_triangle);
    for(int i = 0; i < obs->cf_size; i++) {
        ASSERT(obs->condensed_form[i] == exp->condensed_form[i]);
    }
    for(int i = 0; i < obs->n_samples; i++)
        ASSERT(strcmp(obs->sample_ids[i], exp->sample_ids[i]) == 0);
 
    // error checking
    pm1->stripe_start = 0;
    pm1->stripe_stop = 3;
    pm1->stripe_total = 9;
    pm1->is_upper_triangle = true;

    pm2->stripe_start = 3;
    pm2->stripe_stop = 5;
    pm2->stripe_total = 9;
    pm2->is_upper_triangle = true;
    
    partial_mat_t* pm3 = (partial_mat_t*)malloc(sizeof(partial_mat_t));
    pm3 = make_test_pm();
    pm3->stripe_start = 6;
    pm3->stripe_stop = 9;
    pm3->stripe_total = 9;
    pm3->is_upper_triangle = true;

    exp = mat_three_rep();

    pms[2] = pm1;
    pms[0] = pm2;
    pms[1] = pm3;

    err = merge_partial(pms, 3, 1, &obs);
    ASSERT(err == incomplete_stripe_set);

    pm2->stripe_start = 2;
    pm2->stripe_stop = 6;
    err = merge_partial(pms, 3, 1, &obs);
    ASSERT(err == stripes_overlap);

    pm2->stripe_start = 3;
    pm2->sample_ids[2][0] = 'X';
    err = merge_partial(pms, 3, 1, &obs);
    ASSERT(err == sample_id_consistency);

    pm2->sample_ids[2][0] = 'C';
    pm3->n_samples = 2;
    err = merge_partial(pms, 3, 1, &obs);
    ASSERT(err == partials_mismatch);

    pm3->n_samples = 6;
    pm3->stripe_total = 12;
    err = merge_partial(pms, 3, 1, &obs);
    ASSERT(err == partials_mismatch);
    
    pm3->is_upper_triangle = false;
    pm3->stripe_total = 9;
    err = merge_partial(pms, 3, 1, &obs);
    ASSERT(err == square_mismatch);
    
    SUITE_END();
}

int main(int argc, char** argv) {
    /* one_off and partial are executed as integration tests */    

    //test_write_mat();
    //test_read_mat();
    test_read_write_partial_mat();
    test_merge_partial_mat();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);
  
    printf("\n THE END.\n");
    
    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

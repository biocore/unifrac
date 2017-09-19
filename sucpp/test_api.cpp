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
void test_read_write_partial_mat() {
    SUITE_START("test read/write partial_mat_t");

    partial_mat_t pm;
    pm.n_samples = 5;
    pm.sample_ids = (char**)malloc(sizeof(char*) * 5);
    pm.sample_ids[0] = (char*)malloc(sizeof(char) * 2);
    pm.sample_ids[0][0] = 'A'; pm.sample_ids[0][1] = '\0';
    pm.sample_ids[1] = (char*)malloc(sizeof(char) * 2);
    pm.sample_ids[1][0] = 'B'; pm.sample_ids[1][1] = '\0';
    pm.sample_ids[2] = (char*)malloc(sizeof(char) * 3);
    pm.sample_ids[2][0] = 'C'; pm.sample_ids[2][1] = 'x'; pm.sample_ids[2][2] = '\0';
    pm.sample_ids[3] = (char*)malloc(sizeof(char) * 2);
    pm.sample_ids[3][0] = 'D'; pm.sample_ids[3][1] = '\0';
    pm.sample_ids[4] = (char*)malloc(sizeof(char) * 2);
    pm.sample_ids[4][0] = 'E'; pm.sample_ids[4][1] = '\0';

    pm.stripes = (double**)malloc(sizeof(double*) * 3);
    pm.stripes[0] = (double*)malloc(sizeof(double) * 5);
    pm.stripes[0][0] = 1; pm.stripes[0][1] = 2; pm.stripes[0][2] = 3; pm.stripes[0][3] = 4; pm.stripes[0][4] = 5;
    pm.stripes[1] = (double*)malloc(sizeof(double) * 5);
    pm.stripes[1][0] = 6; pm.stripes[1][1] = 7; pm.stripes[1][2] = 8; pm.stripes[1][3] = 9; pm.stripes[1][4] = 10;
    pm.stripes[2] = (double*)malloc(sizeof(double) * 5);
    pm.stripes[2][0] = 11; pm.stripes[2][1] = 12; pm.stripes[2][2] = 13; pm.stripes[2][3] = 14; pm.stripes[2][4] = 15;
    pm.stripe_start = 0;
    pm.stripe_stop = 3;
    pm.stripe_total = 3;
    pm.is_upper_triangle = true;
    
    io_status err = write_partial("/tmp/ssu_io.dat", &pm);
    ASSERT(err == write_okay);

    partial_mat_t *obs = NULL;
    err = read_partial("/tmp/ssu_io.dat", &obs);
    
    ASSERT(err == read_okay);
    ASSERT(obs->n_samples == 5);
    ASSERT(obs->stripe_start == 0);
    ASSERT(obs->stripe_stop == 3);
    ASSERT(obs->stripe_total == 3);
    ASSERT(strcmp(obs->sample_ids[0], "A") == 0);
    ASSERT(strcmp(obs->sample_ids[1], "B") == 0);
    ASSERT(strcmp(obs->sample_ids[2], "Cx") == 0);
    ASSERT(strcmp(obs->sample_ids[3], "D") == 0);
    ASSERT(strcmp(obs->sample_ids[4], "E") == 0);

    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 5; j++) {
            ASSERT(obs->stripes[i][j] == ((i * 5) + j + 1));
        }
    }
    SUITE_END();
}

//void test_merge_partial_mat() {
//    SUITE_START("test merge partial_mat_t");
//    SUITE_END();
//}

int main(int argc, char** argv) {
    /* one_off and partial are executed as integration tests */    

    //test_write_mat();
    //test_read_mat();
    test_read_write_partial_mat();
    //test_merge_partial_mat();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);
  
    printf("\n THE END.\n");
    
    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

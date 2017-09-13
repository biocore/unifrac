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
void test_write_partial_mat() {
    SUITE_START("test write partial_mat_t");
    SUITE_END();
}

void test_read_partial_mat() {
    SUITE_START("test read partial_mat_t");
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
    test_write_partial_mat();
    test_read_partial_mat();
    //test_merge_partial_mat();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);
  
    printf("\n THE END.\n");
    
    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

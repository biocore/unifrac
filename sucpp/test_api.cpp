#include <iostream>
#include "api.hpp"
#include <cmath>
#include <unordered_set>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <H5Cpp.h>
#include <H5Dpublic.h>

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

template<class TMat, class TReal>
void fill_test_pm(TMat* pm, int case_id) {
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

    if (case_id==0) {
      pm->stripe_start = 0;
      pm->stripe_stop = 3;
      pm->stripe_total = 3;
      pm->stripes = (TReal**)malloc(sizeof(TReal*) * 3);
      pm->stripes[0] = (TReal*)malloc(sizeof(TReal) * 6);
      pm->stripes[0][0] = 1; pm->stripes[0][1] = 2; pm->stripes[0][2] = 3; pm->stripes[0][3] = 4; pm->stripes[0][4] = 5; pm->stripes[0][5] = 6;
      pm->stripes[1] = (TReal*)malloc(sizeof(TReal) * 6);
      pm->stripes[1][0] = 7; pm->stripes[1][1] = 8; pm->stripes[1][2] = 9; pm->stripes[1][3] = 10; pm->stripes[1][4] = 11; pm->stripes[1][5] = 12;
      pm->stripes[2] = (TReal*)malloc(sizeof(TReal) * 6);
      pm->stripes[2][0] = 13; pm->stripes[2][1] = 14; pm->stripes[2][2] = 15; pm->stripes[2][3] = 16; pm->stripes[2][4] = 17; pm->stripes[2][5] = 18;
    } else if (case_id==1) {
      pm->stripe_start = 0;
      pm->stripe_stop = 2;
      pm->stripe_total = 3;
      pm->stripes = (TReal**)malloc(sizeof(TReal*) * 2);
      pm->stripes[0] = (TReal*)malloc(sizeof(TReal) * 6);
      pm->stripes[0][0] = 1; pm->stripes[0][1] = 2; pm->stripes[0][2] = 3; pm->stripes[0][3] = 4; pm->stripes[0][4] = 5; pm->stripes[0][5] = 6;
      pm->stripes[1] = (TReal*)malloc(sizeof(TReal) * 6);
      pm->stripes[1][0] = 7; pm->stripes[1][1] = 8; pm->stripes[1][2] = 9; pm->stripes[1][3] = 10; pm->stripes[1][4] = 11; pm->stripes[1][5] = 12;
    } else { // assume 2
      pm->stripe_start = 2;
      pm->stripe_stop = 3;
      pm->stripe_total = 3;
      pm->stripes = (TReal**)malloc(sizeof(TReal*) * 1);
      pm->stripes[0] = (TReal*)malloc(sizeof(TReal) * 6);
      pm->stripes[0][0] = 16; pm->stripes[0][1] = 17; pm->stripes[0][2] = 18; pm->stripes[0][3] = 16; pm->stripes[0][4] = 17; pm->stripes[0][5] = 18;
    }
    pm->is_upper_triangle = true;
}

partial_mat_t* make_test_pm(int case_id) {
    partial_mat_t* pm = (partial_mat_t*)malloc(sizeof(partial_mat_t));

    fill_test_pm<partial_mat_t,double>(pm,case_id);
    return pm;
}

partial_dyn_mat_t* make_test_pdm(int case_id) {
    partial_dyn_mat_t* pm = (partial_dyn_mat_t*)malloc(sizeof(partial_dyn_mat_t));
    fill_test_pm<partial_dyn_mat_t,double>(pm,case_id);
    pm->offsets = (uint64_t*)calloc(pm->stripe_stop-pm->stripe_start,sizeof(uint64_t));
    pm->filename = strdup("dummy");

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

template<class TMat, class TReal>
TMat* mat_full_three_rep() {
    TMat* res = (TMat*)malloc(sizeof(TMat));
    res->n_samples = 6;
    res->flags=0;
    res->matrix = (TReal*)malloc(sizeof(TReal) * 36);
    TReal * m=res->matrix ;
    m[ 0] =  0; m[ 1] =  1; m[ 2] =  7; m[ 3] = 16; m[ 4] = 11; m[ 5] =  6;
    m[ 6] =  1; m[ 7] =  0; m[ 8] =  2; m[ 9] =  8; m[10] = 17; m[11] = 12;
    m[12] =  7; m[13] =  2; m[14] =  0; m[15] =  3; m[16] =  9; m[17] = 18;
    m[18] = 16; m[19] =  8; m[20] =  3; m[21] =  0; m[22] =  4; m[23] = 10;
    m[24] = 11; m[25] = 17; m[26] =  9; m[27] =  4; m[28] =  0; m[29] =  5;
    m[30] =  6; m[31] = 12; m[32] = 18; m[33] = 10; m[34] =  5; m[35] =  0;

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

    partial_mat_t* pm = make_test_pm(0);
    
    io_status err = write_partial("/tmp/ssu_io.dat", pm);
    ASSERT(err == write_okay);

    {
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

      destroy_partial_mat(&obs);
    }

    {
      partial_dyn_mat_t *obs = NULL;
      err = read_partial_header("/tmp/ssu_io.dat", &obs);
   
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
        ASSERT(obs->stripes[i]==NULL);
      }

      err = read_partial_one_stripe(obs,1);
      ASSERT(err == read_okay);

      ASSERT(obs->stripes[0]==NULL);
      ASSERT(obs->stripes[1]!=NULL);
      ASSERT(obs->stripes[2]==NULL);

      {
        const int i = 1;
        for(int j = 0; j < 6; j++) {
            ASSERT(obs->stripes[i][j] == ((i * 6) + j + 1));
        }
      }

      err = read_partial_one_stripe(obs,0);
      ASSERT(err == read_okay);

      ASSERT(obs->stripes[0]!=NULL);
      ASSERT(obs->stripes[1]!=NULL);
      ASSERT(obs->stripes[2]==NULL);

      {
        const int i = 0;
        for(int j = 0; j < 6; j++) {
            ASSERT(obs->stripes[i][j] == ((i * 6) + j + 1));
        }
      }
    
      err = read_partial_one_stripe(obs,2);
      ASSERT(err == read_okay);


      for(int i = 0; i < 3; i++) {
        ASSERT(obs->stripes[i]!=NULL);
        for(int j = 0; j < 6; j++) {
            ASSERT(obs->stripes[i][j] == ((i * 6) + j + 1));
        }
      }

      destroy_partial_dyn_mat(&obs);
    }

    unlink("/tmp/ssu_io.dat");

    SUITE_END();
}

void test_merge_partial_mat() {
    SUITE_START("test merge partial_mat_t");

    // the easy test
    partial_mat_t* pm1 = make_test_pm(1);
    partial_mat_t* pm2 = make_test_pm(2);

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
    for(unsigned int i = 0; i < obs->cf_size; i++) {
        ASSERT(obs->condensed_form[i] == exp->condensed_form[i]);
    }
    for(unsigned int i = 0; i < obs->n_samples; i++)
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
    for(unsigned int i = 0; i < obs->cf_size; i++) {
        ASSERT(obs->condensed_form[i] == exp->condensed_form[i]);
    }
    for(unsigned int i = 0; i < obs->n_samples; i++)
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
    
    partial_mat_t* pm3 = make_test_pm(0);
    pm3->stripe_start = 6;
    pm3->stripe_stop = 9;
    pm3->stripe_total = 9;

    partial_mat_t* pms_err[3];

    pms_err[2] = pm1;
    pms_err[0] = pm2;
    pms_err[1] = pm3;

    err = merge_partial(pms_err, 3, 1, &obs);
    ASSERT(err == incomplete_stripe_set);

    pm2->stripe_start = 2;
    pm2->stripe_stop = 6;
    err = merge_partial(pms_err, 3, 1, &obs);
    ASSERT(err == stripes_overlap);

    pm2->stripe_start = 3;
    pm2->sample_ids[2][0] = 'X';
    err = merge_partial(pms_err, 3, 1, &obs);
    ASSERT(err == sample_id_consistency);

    pm2->sample_ids[2][0] = 'C';
    pm3->n_samples = 2;
    err = merge_partial(pms_err, 3, 1, &obs);
    ASSERT(err == partials_mismatch);

    pm3->n_samples = 6;
    pm3->stripe_total = 12;
    err = merge_partial(pms_err, 3, 1, &obs);
    ASSERT(err == partials_mismatch);
    
    pm3->is_upper_triangle = false;
    pm3->stripe_total = 9;
    err = merge_partial(pms_err, 3, 1, &obs);
    ASSERT(err == square_mismatch);
    
    SUITE_END();
}

void test_merge_partial_dyn_mat() {
    SUITE_START("test merge partial_dyn_mat_t");

    // the easy test
    partial_dyn_mat_t* pm1 = make_test_pdm(1);
    partial_dyn_mat_t* pm2 = make_test_pdm(2);

    mat_full_fp64_t* exp = mat_full_three_rep<mat_full_fp64_t,double>();

    partial_dyn_mat_t* pms[2];
    pms[0] = pm1;
    pms[1] = pm2;

    mat_full_fp64_t* obs = NULL;
    merge_status err = merge_partial_to_matrix(pms, 2, &obs);
    ASSERT(err == merge_okay);
    ASSERT(obs->n_samples == exp->n_samples);
    for(unsigned int i = 0; i < (obs->n_samples*obs->n_samples); i++) {
        ASSERT(obs->matrix[i] == exp->matrix[i]);
    }
    for(unsigned int i = 0; i < obs->n_samples; i++) {
        ASSERT(strcmp(obs->sample_ids[i], exp->sample_ids[i]) == 0);
    }
    // out of order test

    // recreate deallocated stripes
    ASSERT(pm1->stripes[0]==NULL);
    ASSERT(pm1->stripes[1]==NULL);
    ASSERT(pm2->stripes[0]==NULL);
    destroy_partial_dyn_mat(&pm1);
    destroy_partial_dyn_mat(&pm2);


    pm1 = make_test_pdm(1);
    pm2 = make_test_pdm(2);

    pms[0] = pm2;
    pms[1] = pm1;

    mat_full_fp32_t* exp2 = mat_full_three_rep<mat_full_fp32_t,float>();
    mat_full_fp32_t *obs2 = NULL;
    err = merge_partial_to_matrix_fp32(pms, 2, &obs2);
    ASSERT(err == merge_okay);
    ASSERT(obs2->n_samples == exp2->n_samples);
    for(unsigned int i = 0; i < (obs2->n_samples*obs2->n_samples); i++) {
        ASSERT(obs2->matrix[i] == exp2->matrix[i]);
    }
    for(unsigned int i = 0; i < obs2->n_samples; i++)
        ASSERT(strcmp(obs2->sample_ids[i], exp2->sample_ids[i]) == 0);

 
    ASSERT(pm2->stripes[0]==NULL);
    ASSERT(pm1->stripes[0]==NULL);
    ASSERT(pm1->stripes[1]==NULL);
    destroy_partial_dyn_mat(&pm1);
    destroy_partial_dyn_mat(&pm2);


    pm1 = make_test_pdm(1);
    pm2 = make_test_pdm(2);


    // error checking
    pm1->stripe_start = 0;
    pm1->stripe_stop = 3;
    pm1->stripe_total = 9;
    pm1->is_upper_triangle = true;

    pm2->stripe_start = 3;
    pm2->stripe_stop = 5;
    pm2->stripe_total = 9;
    pm2->is_upper_triangle = true;
    
    partial_dyn_mat_t* pm3 = make_test_pdm(0);
    pm3->stripe_start = 6;
    pm3->stripe_stop = 9;
    pm3->stripe_total = 9;

    partial_dyn_mat_t* pms_err[3];

    pms_err[2] = pm1;
    pms_err[0] = pm2;
    pms_err[1] = pm3;

    err = merge_partial_to_matrix(pms_err, 3, &obs);
    ASSERT(err == incomplete_stripe_set);

    pm2->stripe_start = 2;
    pm2->stripe_stop = 6;
    err = merge_partial_to_matrix(pms_err, 3, &obs);
    ASSERT(err == stripes_overlap);

    pm2->stripe_start = 3;
    pm2->sample_ids[2][0] = 'X';
    err = merge_partial_to_matrix(pms_err, 3, &obs);
    ASSERT(err == sample_id_consistency);

    pm2->sample_ids[2][0] = 'C';
    pm3->n_samples = 2;
    err = merge_partial_to_matrix(pms_err, 3, &obs);
    ASSERT(err == partials_mismatch);

    pm3->n_samples = 6;
    pm3->stripe_total = 12;
    err = merge_partial_to_matrix(pms_err, 3, &obs);
    ASSERT(err == partials_mismatch);
    
    /*
     * Disable for now... not dealing properly with is_upper_triangle == false

    pm3->is_upper_triangle = false;
    pm3->stripe_total = 9;
    err = merge_partial_to_matrix(pms_err, 3, &obs);
    ASSERT(err == square_mismatch);
    */
   
    destroy_mat_full_fp64(&obs);
    // note, we cannot cleanly destroy the partial_dyn_mat_t structures that have been hacked by hand

    SUITE_END();
}

void test_merge_partial_io() {
    SUITE_START("test merge partial_io");

    // the easy test
    partial_mat_t* s1 = make_test_pm(1);
    partial_mat_t* s2 = make_test_pm(2);

    io_status ierr;

    ierr = write_partial("/tmp/ssu_io_1.dat", s1);
    ASSERT(ierr == write_okay);
 
    ierr = write_partial("/tmp/ssu_io_2.dat", s2);
    ASSERT(ierr == write_okay);

    partial_dyn_mat_t* pm1 = NULL;
    partial_dyn_mat_t* pm2 = NULL;

    ierr = read_partial_header("/tmp/ssu_io_1.dat", &pm1);
    ASSERT(ierr == read_okay);

    ierr = read_partial_header("/tmp/ssu_io_2.dat", &pm2);
    ASSERT(ierr == read_okay);

    mat_full_fp64_t* exp = mat_full_three_rep<mat_full_fp64_t,double>();

    partial_dyn_mat_t* pms[2];
    pms[0] = pm1;
    pms[1] = pm2;

    mat_full_fp64_t* obs = NULL;
    merge_status err = merge_partial_to_matrix(pms, 2, &obs);
    ASSERT(err == merge_okay);
    ASSERT(obs->n_samples == exp->n_samples);
    for(unsigned int i = 0; i < (obs->n_samples*obs->n_samples); i++) {
        ASSERT(obs->matrix[i] == exp->matrix[i]);
    }
    for(unsigned int i = 0; i < obs->n_samples; i++) {
        ASSERT(strcmp(obs->sample_ids[i], exp->sample_ids[i]) == 0);
    }
    ASSERT(pm1->stripes[0]==NULL);
    ASSERT(pm1->stripes[1]==NULL);
    ASSERT(pm2->stripes[0]==NULL);

    destroy_mat_full_fp64(&obs);
    destroy_partial_dyn_mat(&pm1);
    destroy_partial_dyn_mat(&pm2);

    // out of order test
    partial_dyn_mat_t* pm1b = NULL;
    partial_dyn_mat_t* pm2b = NULL;

    ierr = read_partial_header("/tmp/ssu_io_1.dat", &pm1b);
    ASSERT(ierr == read_okay);
    
    ierr = read_partial_header("/tmp/ssu_io_2.dat", &pm2b);
    ASSERT(ierr == read_okay);

    pms[0] = pm2b;
    pms[1] = pm1b;

    mat_full_fp32_t* exp2 = mat_full_three_rep<mat_full_fp32_t,float>();
    mat_full_fp32_t *obs2 = NULL;
    err = merge_partial_to_matrix_fp32(pms, 2, &obs2);
    ASSERT(err == merge_okay);
    ASSERT(obs2->n_samples == exp2->n_samples);
    for(unsigned int i = 0; i < (obs2->n_samples*obs2->n_samples); i++) {
        ASSERT(obs2->matrix[i] == exp2->matrix[i]);
    }
    for(unsigned int i = 0; i < obs2->n_samples; i++)
        ASSERT(strcmp(obs2->sample_ids[i], exp2->sample_ids[i]) == 0);


    ASSERT(pm2b->stripes[0]==NULL);
    ASSERT(pm1b->stripes[0]==NULL);
    ASSERT(pm1b->stripes[1]==NULL);

    destroy_mat_full_fp32(&obs2);
    destroy_partial_dyn_mat(&pm1b);
    destroy_partial_dyn_mat(&pm2b);

    unlink("/tmp/ssu_io_1.dat");
    unlink("/tmp/ssu_io_2.dat");

    SUITE_END();
}

void test_merge_partial_mmap() {
    SUITE_START("test merge partial_mmap");

    // the easy test
    partial_mat_t* s1 = make_test_pm(1);
    partial_mat_t* s2 = make_test_pm(2);

    io_status ierr;

    ierr = write_partial("/tmp/ssu_io_1.dat", s1);
    ASSERT(ierr == write_okay);

    ierr = write_partial("/tmp/ssu_io_2.dat", s2);
    ASSERT(ierr == write_okay);

    partial_dyn_mat_t* pm1 = NULL;
    partial_dyn_mat_t* pm2 = NULL;

    ierr = read_partial_header("/tmp/ssu_io_1.dat", &pm1);
    ASSERT(ierr == read_okay);

    ierr = read_partial_header("/tmp/ssu_io_2.dat", &pm2);
    ASSERT(ierr == read_okay);

    mat_full_fp64_t* exp = mat_full_three_rep<mat_full_fp64_t,double>();

    partial_dyn_mat_t* pms[2];
    pms[0] = pm1;
    pms[1] = pm2;

    mat_full_fp32_t* obs = NULL;
    merge_status err = merge_partial_to_mmap_matrix_fp32(pms, 2, "/tmp", &obs);
    ASSERT(err == merge_okay);
    ASSERT(obs->n_samples == exp->n_samples);
    ASSERT(obs->flags != 0);
    for(unsigned int i = 0; i < (obs->n_samples*obs->n_samples); i++) {
        ASSERT(obs->matrix[i] == exp->matrix[i]);
    }
    for(unsigned int i = 0; i < obs->n_samples; i++) {
        ASSERT(strcmp(obs->sample_ids[i], exp->sample_ids[i]) == 0);
    }
    ASSERT(pm1->stripes[0]==NULL);
    ASSERT(pm1->stripes[1]==NULL);
    ASSERT(pm2->stripes[0]==NULL);

    destroy_mat_full_fp32(&obs);
    destroy_partial_dyn_mat(&pm1);
    destroy_partial_dyn_mat(&pm2);

    // test failure due to FS problems

    ierr = read_partial_header("/tmp/ssu_io_1.dat", &pm1);
    ASSERT(ierr == read_okay);

    ierr = read_partial_header("/tmp/ssu_io_2.dat", &pm2);
    ASSERT(ierr == read_okay);

    pms[0] = pm1;
    pms[1] = pm2;


    err = merge_partial_to_mmap_matrix_fp32(pms, 2, "/santa/goes/skiing", &obs);
    ASSERT(err != merge_okay);
    destroy_partial_dyn_mat(&pm1);
    destroy_partial_dyn_mat(&pm2);

    destroy_partial_mat(&s1);
    destroy_partial_mat(&s2);

    unlink("/tmp/ssu_io_1.dat");
    unlink("/tmp/ssu_io_2.dat");


    SUITE_END();
}

void test_to_file_one(const char *method) {

    static const char h5name[]="/tmp/ssu_t1.h5";
    struct stat sbuf;

    ComputeStatus urc;
    int frc;

    // ensure file does not already exist
    frc=stat(h5name,&sbuf);
    if (frc == 0) {
      unlink(h5name);
      frc=stat(h5name,&sbuf);
    }
    ASSERT(frc != 0);

    urc=unifrac_to_file("test.biom","test.tre",h5name,method,false,1.0,false,1,"hdf5",0,NULL);
    ASSERT(urc == okay);

    // first, we check it does exist
    frc=stat(h5name,&sbuf);
    ASSERT(frc == 0);

    {
      try {
        H5::H5File file(h5name, H5F_ACC_RDONLY);
        H5::DataSet mds(file.openDataSet("matrix"));
        H5::DataSpace dataspace(mds.getSpace());

        ASSERT(dataspace.isSimple() == true);
        ASSERT(dataspace.getSimpleExtentNdims() == 2);

        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims, NULL);
        ASSERT(dims[0] == 6);
        ASSERT(dims[1] == 6);
      } catch(...) {
        int rc=1;
        ASSERT(rc == 0); // if we get here is always an error, just to get a nice message
      }
    }

    unlink(h5name);

}

void test_to_file() {
    SUITE_START("test unifrac_to_file");

    test_to_file_one("unweighted");
    test_to_file_one("unweighted_fp32");
    test_to_file_one("weighted_normalized");
    test_to_file_one("weighted_normalized_fp32");
    test_to_file_one("weighted_unnormalized");
    test_to_file_one("weighted_unnormalized_fp32");
    test_to_file_one("generalized");
    test_to_file_one("generalized_fp32");

    SUITE_END();
}

int main(int argc, char** argv) {
    /* one_off and partial are executed as integration tests */    

    //test_write_mat();
    //test_read_mat();
    test_read_write_partial_mat();
    test_merge_partial_mat();
    test_merge_partial_dyn_mat();
    test_merge_partial_io();
    test_merge_partial_mmap();
    test_to_file();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);
  
    printf("\n THE END.\n");
    
    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

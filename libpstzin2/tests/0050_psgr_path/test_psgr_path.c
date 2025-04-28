#define PROG_NAME "test_gr_path"
#define PROG_DESC "checks the routines in {psgr.h}, {psgr_plot.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-04-27 11:33:45 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_gr_path_C_COPYRIGHT \
  "Copyright © 2010  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  test_encode_gamma(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2024-12-23 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_gr_path_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program checks the functions {psgr_path.h}."

#define PROG_INFO_OPTS \
  ""

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <vec.h>
#include <argparser.h>
#include <epswr.h>
#include <r2.h>
#include <i2.h>

#include <psgr_path.h>

typedef struct tepa_options_t
  { int32_t DUMMY;  /* Placeholder. */
  } tepa_options_t;

tepa_options_t *tepa_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {tepa_options_t} record. */

void tepa_test_all(uint32_t n);

psgr_path_t test__psgr_path_throw(i2_t *ip0, uint32_t n, i2_t *ip1);
void test__psgr_path_reverse(psgr_path_t P);
void test__psgr_path_write__psgr_path_read(psgr_path_t P);
void test__psgr_path_center__psgr_path_central_dir__psgr_path_starting_dir(i2_t *ip0, psgr_path_t P, i2_t *ip1);
void test__psgr_path_copy__psgr_path_equal(psgr_path_t P);
void test__psgr_path_concatenate(uint32_t n, i2_t *ip0, i2_t *ipm, i2_t *ip1);
void test__psgr_path_in_star_wedge(i2_t *ip0, i2_t *ipm, i2_t *ip1);
void test__psgr_path_free(psgr_path_t P);

void test__psgr_path_reverse__psgr_path_write(void);

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    tepa_options_t *o = tepa_parse_options(argc, argv);
    
    tepa_test_all(0);
    tepa_test_all(1);
    tepa_test_all(2);
    tepa_test_all(3);
    tepa_test_all(10);

    free(o); o = NULL;
    return 0;
  }
    
void tepa_test_all(uint32_t n)
  { 
    fprintf(stderr, "=== %s (n = %d) ===\n", __FUNCTION__, n);
    i2_t ip0 = (i2_t){{ 3, 4 }};
    i2_t ip1 = (i2_t){{ 15, 17 }};
    i2_t ipm = (i2_t){{ 5, 19 }};

    psgr_path_t P = test__psgr_path_throw(&ip0, n, &ip1);

    test__psgr_path_reverse(P);
    test__psgr_path_write__psgr_path_read(P);
    test__psgr_path_center__psgr_path_central_dir__psgr_path_starting_dir(&ip0, P, &ip1);
    test__psgr_path_copy__psgr_path_equal(P);
    test__psgr_path_concatenate(n, &ip0, &ipm, &ip1);
    test__psgr_path_in_star_wedge(&ip0, &ipm, &ip1);

    /* Cleanup: */
    test__psgr_path_free(P);
  }
    
psgr_path_t test__psgr_path_throw(i2_t *ip0, uint32_t n, i2_t *ip1)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    psgr_path_t P = psgr_path_throw(ip0, n, ip1);
    
    /* Check the position of nodes: */
    r2_t p0 = (r2_t){{ ip0->c[0], ip0->c[1] }};
    r2_t p1 = (r2_t){{ ip1->c[0], ip1->c[1] }};
    r2_t u; r2_sub(&p1, &p0, &u); 
    double d01 = r2_dir(&u, &u);
    double d_step = d01/(n+1);
    double d_prev = - d_step;
    for (int32_t k = -1; k <= (int32_t)n; k++)
      { r2_t pk;
        if (k == -1)
          { pk = p0; }
        else if (k == n)
          { pk = p1; }
        else
          { pk = P.v[k]; }
        r2_t p0k; r2_sub(&pk, &p0, &p0k);
        double d_cmp = r2_dot(&p0k, &u);
        double d_exp = d_prev + d_step;
        demand(fabs(d_cmp - d_exp) < 1.0e-6, "bad spacing of path nodes");
        d_prev = d_cmp;
      }
    assert(fabs(d_prev - d01) < 1.0e-6);
    return P;
  }
    
void test__psgr_path_reverse(psgr_path_t P)
  { fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    psgr_path_t R = psgr_path_reverse(P);
    demand(P.n == R.n, "path_reverse bug (n)");
    demand(P.reverse == ! R.reverse, "path_reverse bug (reverse)");
    demand(P.v == R.v, "path_reverse bug (v)");
  }
     
void test__psgr_path_write__psgr_path_read(psgr_path_t P)
  { fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    char *fname = jsprintf("out/path-n%03d.txt", P.n);
    FILE *wr = open_write(fname, FALSE);
    psgr_path_write(wr, P);
    fprintf(wr, "\n");
    fclose(wr);
    FILE *rd = open_read(fname, FALSE);
    psgr_path_t R = psgr_path_read(rd, 64, 48);
    fclose(rd);
    demand(P.n == R.n, "path write/read bug (n)");
    demand(P.reverse == R.reverse, "path write/read bug (reverse)");
    for (uint32_t k = 0; k < P.n; k++)
      { double d = r2_dist(&(P.v[k]), &(R.v[k]));
        demand(d < 0.001, "path write/read bug (v[k] pixel indices)");
      }
    psgr_path_free(R);
  }
    
void test__psgr_path_center__psgr_path_central_dir__psgr_path_starting_dir(i2_t *ip0, psgr_path_t P, i2_t *ip1)
  { fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    r2_t ctr[2];
    r2_t cdir[2];
    r2_t sdir[2];
    
    i2_t iq0 = *ip0;
    i2_t iq1 = *ip1;
    psgr_path_t Q = P;
    for (int32_t s = 0; s <= 1; s++)
      { r2_t q0 = (r2_t){{ iq0.c[0], iq0.c[1] }};
        r2_t q1 = (r2_t){{ iq1.c[0], iq1.c[1] }};
        r2_t u01; r2_sub(&q1, &q0, &u01); (void)r2_dir(&u01, &u01);
        
        ctr[s] = psgr_path_center(&iq0, Q, &iq1);
        if (P.n == 0)
          { /* Center should be the midpoint: */
            r2_t ctr_exp; r2_mix(0.5, &q0, 0.5, &q1, &ctr_exp);
            double d = r2_dist(&(ctr[s]), &(ctr_exp));
            demand(d <= 1.0e-7, "center bug (n = 0)");
          }
        
        cdir[s] = psgr_path_central_dir(&iq0, Q, &iq1);
        demand(fabs(r2_norm(&(cdir[s])) - 1.0) <= 1.0e-7, "central_dir not normalized");
        if (P.n == 0)
          { /* Central dir should be the segment {q0--q1}: */
            double d = r2_dist(&(cdir[s]), &u01);
            if (d > 1.0e-7) 
              { fprintf(stderr, "iq0 = (%d,%d) iq1 = (%d,%d)\n", iq0.c[0], iq0.c[1], iq1.c[0], iq1.c[1]);
                fprintf(stderr, "cdir[s] = (%+11.8f,%+11.8f)\n", cdir[s].c[0], cdir[s].c[1]);
                fprintf(stderr, "u01     = (%+11.8f,%+11.8f)\n", u01.c[0], u01.c[1]);
              }
            demand(d <= 1.0e-7, "central_dir bug (n = 0)");
          }
        
        sdir[s] = psgr_path_starting_dir(&iq0, Q, &iq1);
        demand(fabs(r2_norm(&(sdir[s])) - 1.0) <= 1.0e-7, "starting_dir not normalized");
        if (P.n == 0)
          { /* Starting dir should be the segment {q0--q1}: */
            double d = r2_dist(&(sdir[s]), &u01);
            demand(d <= 1.0e-7, "starting_dir bug (n = 0)");
          }
        else
          { /* Starting dir should be the segment {q0--P.v[0|n-1]}: */
            r2_t v0 = (s == 0 ? P.v[0] : P.v[P.n-1]);
            r2_t sdir_exp; r2_sub(&v0, &q0, &sdir_exp);
            (void)r2_dir(&sdir_exp, &sdir_exp);
            double d = r2_dist(&(sdir[s]), &sdir_exp);
            demand(d <= 1.0e-7, "starting_dir bug (n >= 1)");
          }

        Q = psgr_path_reverse(Q);
        i2_t iqt = iq0; iq0 = iq1; iq1 = iqt;
      }
    /* The center should be invariant under reversal: */
    demand(r2_dist(&(ctr[0]), &(ctr[1])) <= 1.0e-7, "center/reverse bug");
    
    /* The central dir should reverse under reversal: */
    r2_t cdd; r2_add(&(cdir[0]), &(cdir[1]), &cdd);
    demand(r2_norm(&cdd) <= 1.0e-7, "central_dir/reverse bug");
  }

void test__psgr_path_copy__psgr_path_equal(psgr_path_t P)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    psgr_path_t Q = psgr_path_copy(P);
    demand(psgr_path_equal(P, Q), "copy/equal bug");
    
    psgr_path_t Pr = psgr_path_reverse(P);
    psgr_path_t Qr = psgr_path_copy(Pr);
    demand(psgr_path_equal(Pr, Qr), "copy/equal/reverse bug");

    assert(Pr.v == P.v);
    psgr_path_free(Q);
    psgr_path_free(Qr);
  }

void test__psgr_path_concatenate(uint32_t n, i2_t *ip0, i2_t *ipm, i2_t *ip1)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    for (uint32_t s0 = 0; s0 <= 1; s0++)
      { for (uint32_t s1 = 0; s1 <= 1; s1++)
          { psgr_path_t P0 = psgr_path_throw(ip0, n, ipm);
            if (s0 == 1) { P0 = psgr_path_reverse(P0); }
            psgr_path_t P1 = psgr_path_throw(ipm, n+3, ip1);
            if (s1 == 1) { P1 = psgr_path_reverse(P1); }
            r2_t pm = (r2_t){{ ipm->c[0], ipm->c[1] }};
            psgr_path_t R = psgr_path_concatenate(P0, &pm, P1);
            demand(R.n == P0.n + 1 + P1.n, "concatenate bug, {.n}");
            demand(R.reverse == FALSE, "concatenate bug, {.reverse}");

            for (uint32_t k = 0; k < R.n; k++)
              { r2_t vk_exp;
                if (k < P0.n)
                  { uint32_t k0 = k; vk_exp = (P0.reverse ? P0.v[P0.n-1-k0] : P0.v[k0]); }
                else if (k > P0.n)
                  { uint32_t k1 = k - P0.n - 1; vk_exp = (P1.reverse ? P1.v[P1.n-1-k1] : P1.v[k1]); }
                else
                  { vk_exp = pm; }
                double d = r2_dist(&(R.v[k]), &(vk_exp));
                if (d > 0)
                  { fprintf(stderr, "P0 = "); psgr_path_write(stderr, P0); fprintf(stderr, "\n");
                    fprintf(stderr, "pm = "); r2_print(stderr, &pm); fprintf(stderr, "\n");
                    fprintf(stderr, "P1 = "); psgr_path_write(stderr, P1); fprintf(stderr, "\n");
                    fprintf(stderr, "R =  "); psgr_path_write(stderr, R); fprintf(stderr, "\n");
                    fprintf(stderr, "k = %d", k); r2_print(stderr, &vk_exp); fprintf(stderr, "\n");
                  }
                demand(d == 0, "concatenate bug, {.v[k]}");
              }

            psgr_path_free(P0);
            psgr_path_free(P1);
            psgr_path_free(R);
         }
      }
  }
  
void test__psgr_path_in_star_wedge(i2_t *ip0, i2_t *ipm, i2_t *ip1)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    fprintf(stderr, "!! psgr_path_in_star_wedge NOT TESTED\n");
  }
    
void test__psgr_path_free(psgr_path_t P)
  {
    fprintf(stderr, "--- %s ---\n", __FUNCTION__);
    psgr_path_free(P);
  }

tepa_options_t *tepa_parse_options(int32_t argc, char **argv)
  {
    /* INITIALIZATION: */

    /* Start the command line analyzer {pp}: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);

    /* Process "-help" and "-info" options: */
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    tepa_options_t *o = talloc(1, tepa_options_t);

    /* PARSE KEYWORD ARGUMENTS: */
    
    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }

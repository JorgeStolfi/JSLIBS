/* test program */

#include "foi.h"
#include "foifloat.h"
#include "foimisc.h"
#include "foi1g.h"
#include "foi2z.h"
#include "foi2q.h"
#include "foi2qe.h"
#include "interval.h"
#include "iomisc.h"
#include <math.h>
#include <stdio.h>

/*** INTERNAL PROTOTYPES ***/

int main(void);
void foit_test_interval(FILE *ft);
void foit_test_foi(FILE *ft);
void foit_2d_zeros(void);
void foit_2d_quad(void);
void foit_2d_quad_e(void);
void foit_graph(void);

Float foit_f (Float x, Float y);
Interval foit_fv (Interval x, Interval y);
FOIP foit_ff (FOIP x, FOIP y);
extern char *foit_fname;
extern Interval foit_fxd, foit_fyd;

Float foit_g (Float x);
Interval foit_gv (Interval x);
FOIP foit_gf (FOIP x);
extern char *foit_gname;
extern Interval foit_gxd, foit_gyd;

/*** IMPLEMENTATIONS ***/

int main(void)
  {
    FILE *ft = fopen("foitest.out", "w");

    if (ft == NULL)
      { error ("foitest_main: can't open output file"); }

    foi_init();

    /* foit_test_interval(ft); */
    /* foit_test_foi(ft); */

    fprintf(stderr, "foi_stack_top = %p (%d)\n", foi_top(), (unsigned) foi_top());

    foit_graph();
    /* foit_2d_zeros(); */
    /* foit_2d_quad(); */
    foit_2d_quad_e();
    fclose(ft);
    return (0);
  }

void foit_test_interval(FILE *ft)
  {
    Interval x = {0.5, 2.0};
    Interval y = {2.0, 4.0};
    Interval two = {2.0, 2.0};
    static Interval z;

    fprintf(ft, "\n=== Testing the interval package ===\n\n");

    fprintf(ft, "x =         "); iv_print(ft, x); fputc('\n', ft);
    fprintf(ft, "y =         "); iv_print(ft, x); fputc('\n', ft);
    z = iv_add(x, y);
    fprintf(ft, "x + y =     "); iv_print(ft, z); fputc('\n', ft);
    z = iv_neg(x);
    fprintf(ft, "neg(x) =    "); iv_print(ft, z); fputc('\n', ft);
    z = iv_sub(x, y);
    fprintf(ft, "x - y =     "); iv_print(ft, z); fputc('\n', ft);
    z = iv_mul(x, y);
    fprintf(ft, "x * y =     "); iv_print(ft, z); fputc('\n', ft);
    z = iv_inv(x);
    fprintf(ft, "inv(x) =    "); iv_print(ft, z); fputc('\n', ft);
    z = iv_sqrt(two);
    fprintf(ft, "sqrt(2) =   "); iv_print(ft, z); fputc('\n', ft);
    z = iv_sqrt(x);
    fprintf(ft, "sqrt(x) =   "); iv_print(ft, z); fputc('\n', ft);
    z = iv_mul(z,z);
    fprintf(ft, "sqrt(x)^2 = "); iv_print(ft, z); fputc('\n', ft);
    z = iv_affine(x, 0.1, 1.0, 0.00001);
    fprintf(ft, "affine(x,a,b,g) =  "); iv_print(ft, z); fputc('\n', ft);
  }

void foit_test_foi(FILE *ft)
  {
    Interval xv = { 0.5, 2.0};
    Interval yv = { 2.0, 4.0};
    FOIP x, y, z;

    fprintf(ft, "\n=== Testing the FOI package ===\n\n");

    x = foi_from_interval(xv);
    fprintf(ft, "x =        "); foi_print(ft, x); fputc('\n', ft);

    y = foi_from_interval(yv);
    fprintf(ft, "y =        "); foi_print(ft, y); fputc('\n', ft);

    z = foi_const(2.0, 0.1);
    fprintf(ft, "2+/-0.1 =  "); foi_print(ft, z); fputc('\n', ft);

    z = foi_add(x, y);
    fprintf(ft, "x + y =    "); foi_print(ft, z); fputc('\n', ft);
    z = foi_neg(x);
    fprintf(ft, "neg(x) =   "); foi_print(ft, z); fputc('\n', ft);
    z = foi_mul(x, y);
    fprintf(ft, "x * y =    "); foi_print(ft, z); fputc('\n', ft);

    z = foi_affine(x, 10.0, 5.0, 0.1);
    fprintf(ft, "10x+5+/-0.1 =    "); foi_print(ft, z); fputc('\n', ft);

    z = foi_inv(x);
    fprintf(ft, "1/x    =   "); foi_print(ft, z); fputc('\n', ft);
    z = foi_inv(y);
    fprintf(ft, "1/y    =   "); foi_print(ft, z); fputc('\n', ft);
    z = foi_inv(z);
    fprintf(ft, "1/(1/y) =  "); foi_print(ft, z); fputc('\n', ft);
    z = foi_mul(y, foi_inv(y));
    fprintf(ft, "y*(1/y) =  "); foi_print(ft, z); fputc('\n', ft);

    z = foi_sqrt(x);
    fprintf(ft, "sqrt(x) =  "); foi_print(ft, z); fputc('\n', ft);
    z = foi_mul(z, z);
    fprintf(ft, "sqrt(x)^2 =      "); foi_print(ft, z); fputc('\n', ft);

    /*
    z = foi_sub(x, y);
    fprintf(ft, "x - y =    "); foi_print(ft, z); fputc('\n', ft);
    */
  }

/* Test: plot of zeros of a 2d function */

#include "f1.c"

void foit_2d_zeros(void)
  {
    foi2z_plots(
      "foi2z.ps",
      foit_fname,
      foit_f,
      foit_fv,
      foit_ff,
      foit_fxd, foit_fyd,
      foit_fn, 64
    );
  }

void foit_2d_quad(void)
  {
    foi2q_plots(
      "foi2q.ps",
      foit_fname,
      foit_f,
      foit_fv,
      foit_ff,
      foit_fxd, foit_fyd,
      foit_fn, 64
    );
  }

void foit_2d_quad_e(void)
  {
    foi2qe_plots(
      "foi2qe-1.ps", "foi2qe-2.ps", 
      foit_fname,
      foit_f,
      foit_fv,
      foit_ff,
      foit_fxd, foit_fyd,
      foit_fn, 64
    );
  }

/* Test: plot of graph of a function */

#include "g9.c"

void foit_graph(void)
  {
    foi1g_plots(
      "foi1g.ps",
      foit_gname,
      foit_g,
      foit_gv,
      foit_gf,
      foit_gxd, foit_gyd,
      foit_gn, 128
    );
  }


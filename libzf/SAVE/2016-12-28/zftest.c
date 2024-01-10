/* Quick test of zf.h */
/* Last edited on 2016-12-26 22:32:53 by stolfilocal */

#include <affirm.h>
#include <ia.h>
#include <ia_butfly.h>
#include <flt.h>
#include <zf.h>
#include <stdio.h>

typedef struct eval_data_t
  { int num_evals;
    int num_roots;
    int num_intvs;
  } eval_data_t;

/* PROTOTYPES OF TEST FUNCTIONS */

/* These test functions use the interval-slope method
  to compute their output trapezoid {a}. */

void f_affine        (Interval x, Interval *y, ia_butfly_t *a);
  /* A first-degree function: {y = x - 1/3}. */
    
void f_quad_isolated (Interval x, Interval *y, ia_butfly_t *a);
  /* A quadratic function with two simple roots:
    {y = (x - 1/3)(x - 2/3) = (x - 1/2)^2 - 1/36}. */
    
void f_quad_double   (Interval x, Interval *y, ia_butfly_t *a);
  /* A quadratic function with one double root:
    {y = (x - 1/3)^2}. */

/*** INTERNAL PROTOTYPES ***/

int main (void);

void zt_print_kind(FILE *f, zf_kind_t tv);

void test_zf (
    char *title,
    void (*func) (Interval x, Interval *y, ia_butfly_t *a),
    Interval xd,
    Float epsilon,
    Float delta
  );

/*** IMPLEMENTATIONS ***/

#define NTIMES 10000

int main (void)
  { Interval Unit = (Interval) {Zero, One};
    Float Epsilon = 1.0e-6f;
    Float Delta = 1.0e-6f;
    
    ia_init();
    
    ROUND_NEAR;
    fprintf(stderr, "epsilon = "); flt_print(stderr, Epsilon); fprintf(stderr, "\n");
    fprintf(stderr, "delta   = "); flt_print(stderr, Delta);   fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    test_zf("f_affine",        f_affine,        Unit, Epsilon, Delta);
    test_zf("f_quad_isolated", f_quad_isolated, Unit, Epsilon, Delta);
    test_zf("f_quad_double",   f_quad_double,   Unit, Epsilon, Delta);
    
    return(0);
  }
    
void test_zf (
    char *title,
    void (*func) (Interval x, Interval *y, ia_butfly_t *a),
    Interval xd,
    Float epsilon,
    Float delta
  )
  { eval_data_t ed;
    zf_kind_t tn;
    
    fprintf(stderr, "---------------------------------------------\n");
    fprintf(stderr, "testing zf\n");
    fprintf(stderr, "function = %s\n", title);
    fprintf(stderr, "interval = "); ia_print(stderr, xd); fprintf(stderr, "\n");
    fprintf(stderr, "epsilon = "); flt_print(stderr, epsilon); fprintf(stderr, "\n");
    fprintf(stderr, "delta =   "); flt_print(stderr, delta);   fprintf(stderr, "\n");
    
    ed.num_evals = 0;
    ed.num_roots = 0;
    ed.num_intvs = 0;
    
    auto void do_eval(Interval *xr, Interval *yr, ia_butfly_t *a);
      /* The {eval} parameter of {zf_enum_zeros}: 
        evaluates the function {func}, and prints a trace of the call. */

    auto bool_t do_report(Interval *xr, Interval *yr, zf_kind_t kind);
      /* The {report} parameter of the the zero finder: 
         Prints the output interval. */
    
    tn = zf_enum_zeros(
      do_eval, 
      do_report,
      xd,
      epsilon,
      delta
    );

    fprintf(stderr, "num evals =      %d\n", ed.num_evals);
    fprintf(stderr, "intvs reported = %d\n", ed.num_intvs);
    fprintf(stderr, "roots found =    %d\n", ed.num_roots);

    fprintf(stderr, "next interval =  "); zt_print_kind(stderr, tn); fprintf(stderr, "\n");

    return;
    
    auto void do_eval(Interval *xr, Interval *yr, ia_butfly_t *a)
      {  (ed.num_evals)++;
        fprintf(stderr, "  x   = "); ia_print(stderr, *xr);

        /* Evaluate {F} and package the result as a butterfly: */
        func(*xr, yr, a);

        fprintf(stderr, "  f(x) = "); ia_print(stderr, *yr); fprintf(stderr, "\n");
        fprintf(stderr, "  butterfly:\n    ");
        ia_butfly_print(stderr, a, "\n    ");
        fprintf(stderr, "\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "\n");
      }

    bool_t do_report(Interval *xr, Interval *yr, zf_kind_t kind)
      { (ed.num_intvs)++;
        if (kind == zf_kind_root) (ed.num_roots)++;

        fprintf(stderr, "  OUTPUT: kind = ");
        zt_print_kind(stderr, kind);
        fprintf(stderr, "\n");
        fprintf(stderr, "  x   = "); ia_print(stderr, *xr); 
        fprintf(stderr, "  y   = "); ia_print(stderr, *yr);
        fprintf(stderr, "\n");
        fprintf(stderr, "\n");
        return (0);
      }
  }

void zt_print_kind(FILE *f, zf_kind_t kind)
  { switch (kind) 
      { case zf_kind_root: 
	  fprintf(f, "root       ");
          break;
	case zf_kind_positive:
	  fprintf(f, "positive   ");
          break;
	case zf_kind_negative:
	  fprintf(f, "negative   ");
          break;
	case zf_kind_undefined:
	  fprintf(f, "undefined  ");
          break;
	case zf_kind_mixed:
	  fatalerror("zftest: reported a mixed interval");
          break;
	default:
	  fatalerror("zftest: unknown interval kind");
      }
  }

void f_affine (Interval x, Interval *y, ia_butfly_t *a)
  { 
    Interval M1o3 = ia_inv(ia_int_const(-3));
    *y = ia_add(x, M1o3);
    /* Return a butterfly with a trivial lo trapezoid: */
    a->tp[0].x.lo = x.lo;
    a->tp[0].x.hi = a->tp[1].x.lo = x.lo;
    a->tp[1].x.hi = x.hi;
    
    a->tp[0].yxlo = ia_shift(M1o3, x.lo);
    a->tp[0].yxhi = a->tp[1].yxlo = a->tp[0].yxlo;
    a->tp[1].yxhi = ia_shift(M1o3, x.hi);
  }

void f_quad_isolated (Interval x, Interval *y, ia_butfly_t *a)
  { Interval M1o2 = ia_const(-Half, Zero);
    Interval M1o36 = ia_inv(ia_int_const(-36));

    Interval M1o3 = ia_inv(ia_int_const(-3));
    Interval M2o3 = ia_div(ia_int_const(2), ia_int_const(-3));

    *y = ia_add(ia_sqr(ia_add(x, M1o2)), M1o36);

    /* Return a butterfly with two trapezoids: */
    
    auto void eval_trapezoid(Interval x, Interval *yxlo, Interval *yxhi);
      /* Given {x}, sets {yxlo, yxhi} so that the trapezoid {(x,yxlo,yxhi)} 
        encloses the function's graph. */

    a->tp[0].x.lo = x.lo;
    a->tp[0].x.hi = a->tp[1].x.lo = ia_mid(x);
    a->tp[1].x.hi = x.hi;

    eval_trapezoid(a->tp[0].x, &(a->tp[0].yxlo), &(a->tp[0].yxhi));
    eval_trapezoid(a->tp[1].x, &(a->tp[1].yxlo), &(a->tp[1].yxhi));
    
    return;
      
    void eval_trapezoid(Interval xt, Interval *yxlo, Interval *yxhi)
      { /* Start with intervals that contain the graph's endpoints: */
        (*yxlo) = ia_mul(ia_shift(M1o3, xt.lo), ia_shift(M2o3, xt.lo));
        (*yxhi) = ia_mul(ia_shift(M1o3, xt.hi), ia_shift(M2o3, xt.hi));
        /* Now expand the intervals down to enclose the whole graph: */

        ROUND_UP; 
        Float dy = (xt.hi - xt.lo) * (xt.hi - xt.lo) / Four;

        ROUND_DOWN;
        yxlo->lo = yxlo->lo - dy;
        yxhi->lo = yxhi->lo - dy;
      }
  }

void f_quad_double (Interval x, Interval *y, ia_butfly_t *a)
  { Interval M1o3 = ia_inv(ia_int_const(-3));

    *y = ia_sqr(ia_add(x, M1o3));
    
    /* Return a butterfly with two trapezoids: */
    
    auto void eval_trapezoid(Interval x, Interval *yxlo, Interval *yxhi);
      /* Given {x}, sets {yxlo, yxhi} so that the trapezoid {(x,yxlo,yxhi)} 
        encloses the function's graph. */

    a->tp[0].x.lo = x.lo;
    a->tp[0].x.hi = a->tp[1].x.lo = ia_mid(x);
    a->tp[1].x.hi = x.hi;

    eval_trapezoid(a->tp[0].x, &(a->tp[0].yxlo), &(a->tp[0].yxhi));
    eval_trapezoid(a->tp[1].x, &(a->tp[1].yxlo), &(a->tp[1].yxhi));
    
    return;
      
    void eval_trapezoid(Interval xt, Interval *yxlo, Interval *yxhi)
      { /* Start with intervals that contain the graph's endpoints: */
        (*yxlo) = ia_sqr(ia_shift(M1o3, xt.lo));
        (*yxhi) = ia_sqr(ia_shift(M1o3, xt.hi));
        /* Now expand the intervals down to enclose the whole graph: */

        ROUND_UP; 
        Float dy = (xt.hi - xt.lo) * (xt.hi - xt.lo) / Four;

        ROUND_DOWN;
        yxlo->lo = yxlo->lo - dy;
        yxhi->lo = yxhi->lo - dy;
      }
  }

#include <js.h>
#include <ia.h>
#include <flt.h>
#include <zerofind.h>
#include <ioprotos.h>
#include <stdio.h>

typedef struct {
    IntervalPair (*func) (Interval x);
    int num_evals;
    int num_roots;
    int num_intvs;
  } eval_data_t;

/*** INTERNAL PROTOTYPES ***/

int main (void);

IntervalPair f_affine      (Interval x);
IntervalPair f_quad_single (Interval x);
IntervalPair f_quad_double (Interval x);

IntervalPair zerotest_eval (Interval *xv, void *data);
int zerotest_report (Interval *xv, Interval *yv, zf_type tv, void *data);

void zt_print_type(FILE *f, zf_type tv);

void test_zerofind (
    char *title,
    IntervalPair (*func) (Interval x),
    Interval xd,
    Float epsilon,
    Float delta
  );

/*** IMPLEMENTATIONS ***/

#define NTIMES 10000

int main (void)
  {
    Interval Unit = (Interval) {Zero, One};
    Float Epsilon = 1.0e-6;
    Float Delta = 1.0e-6;
    
    ia_init();
    
    ROUND_NEAR;
    fprintf(stderr, "epsilon = "); flt_print(stderr, Epsilon); fprintf(stderr, "\n");
    fprintf(stderr, "delta   = "); flt_print(stderr, Delta);   fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    test_zerofind("f_affine",      f_affine,      Unit, Epsilon, Delta);
    test_zerofind("f_quad_single", f_quad_single, Unit, Epsilon, Delta);
    test_zerofind("f_quad_double", f_quad_double, Unit, Epsilon, Delta);
    
    return(0);
  }
    
void test_zerofind (
    char *title,
    IntervalPair (*func) (Interval x),
    Interval xd,
    Float epsilon,
    Float delta
  )
  {
    eval_data_t ed;
    zf_type tn;
    
    fprintf(stderr, "---------------------------------------------\n");
    fprintf(stderr, "testing zerofind\n");
    fprintf(stderr, "function = %s\n", title);
    fprintf(stderr, "interval = "); ia_print(stderr, xd); fprintf(stderr, "\n");
    fprintf(stderr, "epsilon = "); flt_print(stderr, epsilon); fprintf(stderr, "\n");
    fprintf(stderr, "delta =   "); flt_print(stderr, delta);   fprintf(stderr, "\n");
    
    ed.func = func;
    ed.num_evals = 0;
    ed.num_roots = 0;
    ed.num_intvs = 0;
    
    tn = zerofind(
      zerotest_eval, 
      zerotest_report,
      &ed, 
      xd,
      epsilon,
      delta
    );

    fprintf(stderr, "num evals =      %d\n", ed.num_evals);
    fprintf(stderr, "intvs reported = %d\n", ed.num_intvs);
    fprintf(stderr, "roots found =    %d\n", ed.num_roots);

    fprintf(stderr, "next interval =  "); zt_print_type(stderr, tn); fprintf(stderr, "\n");
  }

IntervalPair zerotest_eval (Interval *xv, void *data)
  { eval_data_t *ed = (eval_data_t *) data;
    IntervalPair y;
    (ed->num_evals)++;
    fprintf(stderr, "  xv = "); ia_print(stderr, *xv);
    y = (ed->func)(*xv);
    fprintf(stderr, " --> a = "); 
    ia_print(stderr, y.a); 
    fprintf(stderr, "  b = "); 
    ia_print(stderr, y.b); 
    fprintf(stderr, "\n");
    return (y);
  }
  
int zerotest_report (Interval *xv, Interval *yv, zf_type tv, void *data)
  { eval_data_t *ed = (eval_data_t *) data;

    (ed->num_intvs)++;
    if (tv == zf_type_root) (ed->num_roots)++;

    fprintf(stderr, "  ** type = ");
    zt_print_type(stderr, tv);
    fprintf(stderr, "  xv = "); 
    ia_print(stderr, *xv); 
    fprintf(stderr, "  yv = ");
    ia_print(stderr, *yv);
    fprintf(stderr, "\n");
    return (0);
  }

void zt_print_type(FILE *f, zf_type tv)
  { switch (tv) 
      {
	case zf_type_root: 
	  fprintf(f, "root       ");
          break;
	case zf_type_positive:
	  fprintf(f, "positive   ");
          break;
	case zf_type_negative:
	  fprintf(f, "negative   ");
          break;
	case zf_type_undefined:
	  fprintf(f, "undefined  ");
          break;
	case zf_type_complex:
	  error("zerotest: reported a complex interval");
          break;
	default:
	  error("zerotest: unknown interval type");
      }
  }

IntervalPair f_affine      (Interval x)
  { /* y = x - 1/3 */
    Interval M1o3 = ia_inv(ia_int_const(-3));
    IntervalPair yp;
    yp.a = ia_shift(M1o3, x.lo);
    yp.b = ia_shift(M1o3, x.hi);
    return (yp);
  }

IntervalPair f_quad_single (Interval x)
  { /* y = (x - 1/3)(x - 2/3) */
    
    Interval M1o3 = ia_inv(ia_int_const(-3));
    Interval M2o3 = ia_div(ia_int_const(2), ia_int_const(-3));
    Float dy;
    IntervalPair yp;
    
    yp.a = ia_mul(ia_shift(M1o3, x.lo), ia_shift(M2o3, x.lo));
    yp.b = ia_mul(ia_shift(M1o3, x.hi), ia_shift(M2o3, x.hi));

    ROUND_UP; dy = (x.hi - x.lo) * (x.hi - x.lo) / Four;

    ROUND_DOWN;
    yp.a.lo = yp.a.lo - dy;
    yp.b.lo = yp.b.lo - dy;
    return (yp);
  }

IntervalPair f_quad_double (Interval x)
  { /* y = (x - 1/3)^2 */
    
    Interval M1o3 = ia_inv(ia_int_const(-3));
    Float dy;
    IntervalPair yp;
    
    yp.a = ia_sqr(ia_shift(M1o3, x.lo));
    yp.b = ia_sqr(ia_shift(M1o3, x.hi));

    ROUND_UP; dy = (x.hi - x.lo) * (x.hi - x.lo) / Four;

    ROUND_DOWN;
    yp.a.lo = yp.a.lo - dy;
    yp.b.lo = yp.b.lo - dy;
    return (yp);
  }

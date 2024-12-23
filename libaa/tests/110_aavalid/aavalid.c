/* Validation of AA ops */
/* Last edited on 2024-12-21 11:22:29 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <flt.h>
#include <ia.h>
#include <affirm.h>
#include <aa.h>
#include <jsrandom.h>

#define DEBUG 0
#define MAXEPS 4

typedef enum {
    aavalid_op_from_interval = 0,
    aavalid_op_int_const =     1,
    aavalid_op_affine =        2,
    aavalid_op_affine_2 =      3,
    aavalid_op_add =     4,
    aavalid_op_sub =     5,
    aavalid_op_neg =     6,
    aavalid_op_inv =     7,
    aavalid_op_mul =     8,
    aavalid_op_div =     9,
    aavalid_op_sqr =    10,
    aavalid_op_sqrt =   11,
    aavalid_op_abs =    12,
    aavalid_op_max =    13,
    aavalid_op_min =    14,
    aavalid_op_cancel1 =  15,
    aavalid_op_cancel2 =  16
  } aavalid_op_t;

typedef struct {  /* Data for one AA operands/result tuple */
    aavalid_op_t  op;  /* Operation code */
    char *op_name;   /* Operation's name */
    /* Miscellaneous arguments: */
    int intx;        /* Argument for {aa_int_const} */
    Interval iax;    /* Argument for {aa_from_interval} */
    /* Arguments for {aa_affine} and {aa_affine_2}: */
    Float alpha;     /* Multiplier for {x} */
    Float beta;      /* Multiplier for {y} (if {aa_affine_2}) */
    Float zeta;      /* Denominator */
    Float gamma;     /* Constant term */
    Float delta;     /* Magnitude of error term */
    /* Main arguments and results: */
    AAP x;           /* First argument of {op} */
    AAP y;           /* Second argument of {op} (if binary) */
    AAP z;           /* Result of {aa_op(x, y)} */
    /* Validation intervals: */
    Interval x_fix;  /* Result of fixing epsilons in {x} */
    Interval y_fix;  /* Result of fixing epsilons in {y} */
    Interval z_fix;  /* Result of fixing epsilons in {z} */
    Interval z_cmp;  /* Result of {ia_op(x_fix, y_fix)} */
  } AAData;

/*** INTERNAL PROTOTYPES ***/

int main(int argc, char *argv[]);

void aavalid_test_op (
    aavalid_op_t op, /* The operation to test */
    char *op_name,   /* The operation's name */
    int nforms,      /* Number of affine evaluations */
    int npoints      /* Number of sample points */
  );
  /* Tests an AA operation {nforms} times.  For each test, generates
     a set {x, y, ...} of AA operands, evaluates {op} on them, and then
     performs {npoints} pointwise checks on the resulting AA form {z}.
     
     Each pointwise check consists of selecting a random assignment of
     values to the noise symbols appearing in the AA operands {x, y,
     ..}; fixing those noise symbols in {x, y, ...}, which results in
     almost point-like intervals {x_fix, y_fix, ...}; applying {op} to those
     intervals, in IA arithmetic, which results in an interval {z_cmp};
     fixing the same noise symbols in the AA result {z}, which results
     in an interval {z_fix}; and finally checking that {z_cmp} and {z_fix}
     have nonempty intersection.  */

void aavalid_throw_and_eval_aa_tuple(AAData *d);
  /* Generates an adequate AA operand tuple {d->x}, {d->y}
    for {d->op}, and performs {aa_op} on it, storing the result
    in {d->z}. */
    
void aavalid_throw_and_fix_eps(AAData *d, AATerm eps[]);
  /* Generates a random assignment of values {eps} for the noise symbols
    {[0..MAXEPS-1]}, and uses it to fix {d->x}, {d->y}, and {d->z},
    obtaining {d->x_fix}, {d->y_fix}, and {d->z_fix}.  Then
    performs {d->op} with IA arithmetic on {d->x_fix} and
    {d->y_fix}, obtaining {d->z_cmp}. */
    
void aavalid_print_data_and_eps(AAData *d, AATerm eps[]);
  /* Prints {*d} and {eps[0..MAXEPS-1]} to {stderr. */
    
float aavalid_throw_eps(void);
  /* Generates a random float in {[-1 .. +1]} with the following distribution:
    0, -1, +1, each with probability 0.25; uniform in {[-1 __ +1]} with 
    probability 0.25. */
    
Float aavalid_throw_float(void);
  /* Generates a random float that is 1, 0, -1, or anything, with
    probability 0.25 each. */

int aavalid_is_zero(AAP x);
  /* TRUE if the AA {x} is exactly zero. */

/*** MAIN PROGRAM ***/

int main(int argc, char *argv[])
  {
    flt_init();
    ia_init();
    aa_init();
  
    aavalid_test_op(aavalid_op_from_interval,   "from_interval", 100,   1);
    aavalid_test_op(aavalid_op_int_const,       "int_const",     100,   1);

    aavalid_test_op(aavalid_op_cancel1,   "cancel1",    200, 5000);
    aavalid_test_op(aavalid_op_cancel2,   "cancel2",    200, 5000);

    aavalid_test_op(aavalid_op_affine,   "affine",    200, 5000);
    aavalid_test_op(aavalid_op_affine_2, "affine_2",  200, 5000);

    aavalid_test_op(aavalid_op_add,      "add",       200, 5000);
    aavalid_test_op(aavalid_op_sub,      "sub",       200, 5000);
    aavalid_test_op(aavalid_op_neg,      "neg",       200, 5000);
    aavalid_test_op(aavalid_op_inv,      "inv",       200, 5000);
    aavalid_test_op(aavalid_op_mul,      "mul",       200, 5000);
    aavalid_test_op(aavalid_op_div,      "div",       200, 5000);
    aavalid_test_op(aavalid_op_sqr,      "sqr",       200, 5000);
    aavalid_test_op(aavalid_op_sqrt,     "sqrt",      200, 5000);
    aavalid_test_op(aavalid_op_abs,      "abs",       200, 5000);
    aavalid_test_op(aavalid_op_max,      "max",       200, 5000);
    aavalid_test_op(aavalid_op_min,      "min",       200, 5000);
    
    return(0);

  }
		     
void aavalid_test_op (
    aavalid_op_t op,   /* The operation to test */
    char *op_name,   /* The operation's name */
    int nforms,      /* Number of affine evalutaions */
    int npoints      /* Number of sample points */
  )
  {
    AAData d;
    AATerm eps[MAXEPS];
    int iform, ipoint;
    int bad, debug;

    MemP frame;
    frame = aa_top();

    fprintf(stderr, "\n");
    fprintf(stderr, "--------------------------------------------------------\n");
    fprintf(stderr, "testing %s\n", op_name);
    fprintf(stderr, "\n");
    
    srand(314159);
    
    /* The following is a trick to ensure that the noise symbols 0 thru {MAXEPS-1},
       which will be used in the test operands, are not used for rounding
       and truncation terms: */
    (void) aa_throw(MAXEPS);
    aa_flush(frame);
    
    d.op = op;
    d.op_name = op_name;
    
    for (iform=0; iform < nforms; iform++)
      {
        aavalid_throw_and_eval_aa_tuple(&d);
        for (ipoint=0; ipoint < npoints; ipoint++)
          { 
            aavalid_throw_and_fix_eps(&d, eps);
            bad = ((d.z_cmp.hi < d.z_fix.lo) || (d.z_cmp.lo > d.z_fix.hi));
            debug = ((ipoint == 0) & (iform < 3));
            if (bad || debug) 
              { fprintf(stderr, "iform = %d  ipoint = %d\n", iform, ipoint);
                aavalid_print_data_and_eps(&d, eps);
                affirm(!bad, "z_fix and z_cmp are disjoint");
                fprintf(stderr, "\n");
              }
          }
	aa_flush(frame);
      }

    fprintf(stderr, "\n");
    fprintf(stderr, "--------------------------------------------------------\n");
  }
  
void aavalid_throw_and_eval_aa_tuple(AAData *d)
  {
    switch (d->op)
      {
	case aavalid_op_from_interval:
	  d->iax = ia_throw();
          d->z = aa_from_interval(d->iax);
	  break;
	case aavalid_op_int_const:
	  d->intx = rand()-rand();
          d->z = aa_int_const(d->intx);
	  break;
	case aavalid_op_affine:
	  d->alpha = aavalid_throw_float();
          do { d->zeta = aavalid_throw_float(); } while (d->zeta == Zero);
          d->gamma = aavalid_throw_float();
          /* d->delta = aavalid_throw_float(); */
          d->delta = Zero;
          d->x = aa_throw(MAXEPS);
          d->y = aa_throw(MAXEPS);
          d->z = aa_affine(d->x, d->alpha, d->zeta, d->gamma, d->delta);
	  break;
	case aavalid_op_affine_2:
	  d->alpha = aavalid_throw_float();
          d->beta = aavalid_throw_float();
          do { d->zeta = aavalid_throw_float(); } while (d->zeta == Zero);
          d->gamma = aavalid_throw_float();
          /* d->delta = aavalid_throw_float(); */
          d->delta = Zero;
          d->x = aa_throw(MAXEPS);
          d->y = aa_throw(MAXEPS);
          d->z = aa_affine_2(d->x, d->alpha, d->y, d->beta, d->zeta, d->gamma, d->delta);
	  break;
	case aavalid_op_add:
	  d->x = aa_throw(MAXEPS);
          d->y = aa_throw(MAXEPS);
          d->z = aa_add(d->x, d->y);
	  break;
	case aavalid_op_sub:
	  d->x = aa_throw(MAXEPS);
          d->y = aa_throw(MAXEPS);
          d->z = aa_sub(d->x, d->y);
	  break;
	case aavalid_op_neg:
	  d->x = aa_throw(MAXEPS);
          d->z = aa_neg(d->x);
	  break;
	case aavalid_op_inv:
	  do { d->x = aa_throw(MAXEPS); } while (aavalid_is_zero(d->x));
          d->z = aa_inv(d->x);
	  break;
	case aavalid_op_mul:
	  d->x = aa_throw(MAXEPS);
          d->y = aa_throw(MAXEPS);
          d->z = aa_mul(d->x, d->y);
	  break;
	case aavalid_op_div:
	  d->x = aa_throw(MAXEPS);
          do { d->y = aa_throw(MAXEPS); } while (aavalid_is_zero(d->y));
          d->z = aa_div(d->x, d->y);
	  break;
	case aavalid_op_sqr:
	  d->x = aa_throw(MAXEPS);
          d->z = aa_sqr(d->x);
	  break;
	case aavalid_op_sqrt:
	  do { d->x = aa_throw(MAXEPS); } while (aa_range(d->x).hi < Zero);
          d->z = aa_sqrt(d->x);
	  break;
	case aavalid_op_abs:
	  d->x = aa_throw(MAXEPS);
          d->z = aa_abs(d->x);
	  break;
	case aavalid_op_max:
	  d->x = aa_throw(MAXEPS);
          d->y = aa_throw(MAXEPS);
          d->z = aa_max(d->x, d->y);
	  break;
	case aavalid_op_min:
	  d->x = aa_throw(MAXEPS);
          d->y = aa_throw(MAXEPS);
          d->z = aa_min(d->x, d->y);
	  break;
	case aavalid_op_cancel1:
	  d->iax = ia_const(-One, flt_random() * flt_random_mag(-23, 2));
          d->y = aa_throw(MAXEPS);
          d->z = aa_add(d->y, aa_mul(d->y, aa_from_interval(d->iax)));
	  break;
	case aavalid_op_cancel2:
	  d->x = aa_throw(MAXEPS);
          { Float a, b;
            ROUND_NEAR; a = (Float)(3.0/7.0); b = (Float)(4.0/7.0);
            d->y = aa_add(aa_affine(d->x, a, One, Zero, Zero), aa_affine(d->x, b, One, Zero, Zero));
          }
          d->z = aa_sub(d->x, d->y);
	  break;
	default:
	  fatalerror("aavalid_throw_and_eval_aa_tuple: bad op code");
      }
  }

void aavalid_throw_and_fix_eps(AAData *d, AATerm eps[])
  {
    int i;
    MemP frame = aa_top();
    
    for (i=0; i<MAXEPS; i++)
      { eps[i].id = i;
        eps[i].coef = aavalid_throw_eps();
      }
      
    switch (d->op)
      {
	case aavalid_op_from_interval:
	  d->z_cmp = d->iax;
	  break;
	case aavalid_op_int_const:
          d->z_cmp = ia_int_const(d->intx);
	  break;
	case aavalid_op_affine:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->z_cmp = ia_affine(d->x_fix, d->alpha, d->zeta, d->gamma, d->delta);
	  break;
	case aavalid_op_affine_2:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_affine_2(d->x_fix, d->alpha, d->y_fix, d->beta, d->zeta, d->gamma, d->delta);
	  break;
	case aavalid_op_add:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_add(d->x_fix, d->y_fix);
	  break;
	case aavalid_op_sub:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_sub(d->x_fix, d->y_fix);
	  break;
	case aavalid_op_neg:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->z_cmp = ia_neg(d->x_fix);
	  break;
	case aavalid_op_inv:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->z_cmp = ia_inv(d->x_fix);
	  break;
	case aavalid_op_mul:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_mul(d->x_fix, d->y_fix);
	  break;
	case aavalid_op_div:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_div(d->x_fix, d->y_fix);
	  break;
	case aavalid_op_sqr:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->z_cmp = ia_sqr(d->x_fix);
	  break;
	case aavalid_op_sqrt:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          if ((d->x_fix).hi < Zero)
            { d->z_cmp = ia_full(); }
          else
            { d->z_cmp = ia_sqrt(d->x_fix); }
	  break;
	case aavalid_op_abs:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->z_cmp = ia_abs(d->x_fix);
	  break;
	case aavalid_op_max:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_max(d->x_fix, d->y_fix);
	  break;
	case aavalid_op_min:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_min(d->x_fix, d->y_fix);
	  break;
	case aavalid_op_cancel1:
	  d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_add(d->y_fix, ia_mul(d->y_fix, d->iax));
	  break;
	case aavalid_op_cancel2:
	  d->x_fix = aa_range(aa_fix_eps(d->x, MAXEPS, eps));
          d->y_fix = aa_range(aa_fix_eps(d->y, MAXEPS, eps));
          d->z_cmp = ia_sub(d->x_fix, d->y_fix);
	  break;
	default:
	  fatalerror("aavalid_throw_and_fix_eps: bad op code");
      }
    d->z_fix = aa_range(aa_fix_eps(d->z, MAXEPS, eps));
    aa_flush(frame);
  }

void aavalid_print_data_and_eps(AAData *d, AATerm eps[])
  {
    int i;
    
    switch (d->op)
      {
	case aavalid_op_from_interval:
	  fprintf(stderr, "iax = "); ia_print(stderr, d->iax); fprintf(stderr, "\n");
	  break;
	case aavalid_op_int_const:
          fprintf(stderr, "intx = %d\n", d->intx);
          { Interval test = ia_int_const(d->intx); 
            fprintf(stderr, "test = "); ia_print(stderr, test); fprintf(stderr, "\n");
          }
	  break;
	case aavalid_op_affine:
	  fprintf(stderr, "alpha = "); flt_print(stderr, d->alpha); fprintf(stderr, "\n");
	  fprintf(stderr, "zeta  = "); flt_print(stderr, d->zeta);  fprintf(stderr, "\n");
	  fprintf(stderr, "gamma = "); flt_print(stderr, d->gamma); fprintf(stderr, "\n");
	  fprintf(stderr, "delta = "); flt_print(stderr, d->delta); fprintf(stderr, "\n");
          fprintf(stderr, "x = "); aa_print(stderr, d->x); fprintf(stderr, "\n");
	  break;
	case aavalid_op_affine_2:
	  fprintf(stderr, "alpha = "); flt_print(stderr, d->alpha); fprintf(stderr, "\n");
	  fprintf(stderr, "beta  = "); flt_print(stderr, d->beta);  fprintf(stderr, "\n");
	  fprintf(stderr, "zeta  = "); flt_print(stderr, d->zeta);  fprintf(stderr, "\n");
	  fprintf(stderr, "gamma = "); flt_print(stderr, d->gamma); fprintf(stderr, "\n");
	  fprintf(stderr, "delta = "); flt_print(stderr, d->delta); fprintf(stderr, "\n");
          fprintf(stderr, "x = "); aa_print(stderr,d->x); fprintf(stderr,"\n");
          fprintf(stderr, "y = "); aa_print(stderr,d->y); fprintf(stderr,"\n");
	  break;
	case aavalid_op_add:
	case aavalid_op_sub:
	case aavalid_op_mul:
	case aavalid_op_div:
	case aavalid_op_max:
	case aavalid_op_min:
        case aavalid_op_cancel2:
          fprintf(stderr, "x = "); aa_print(stderr, d->x); fprintf(stderr, "\n");
          fprintf(stderr, "y = "); aa_print(stderr, d->y); fprintf(stderr, "\n");
	  break;
	case aavalid_op_neg:
	case aavalid_op_inv:
	case aavalid_op_sqr:
	case aavalid_op_sqrt:
	case aavalid_op_abs:
          fprintf(stderr, "x = "); aa_print(stderr,d->x); fprintf(stderr,"\n");
	  break;
	case aavalid_op_cancel1:
	  fprintf(stderr, "iax = "); ia_print(stderr, d->iax); fprintf(stderr, "\n");
          fprintf(stderr, "y = "); aa_print(stderr, d->y); fprintf(stderr, "\n");
	  break;
	default:
	  fatalerror("aavalid_throw_and_fix_eps: bad op code");
      }
    fprintf(stderr,"z = "); aa_print(stderr,d->z); fprintf(stderr,"\n");
    fprintf(stderr, "\n");
    for (i=0; i<MAXEPS; i++)
      {
	fprintf(stderr, "eps[%ld] = ", eps[i].id);
        ROUND_NEAR; flt_print(stderr, eps[i].coef);
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
    switch (d->op)
      {
	case aavalid_op_from_interval:
	case aavalid_op_int_const:
	  break;
	case aavalid_op_affine_2:
	case aavalid_op_add:
	case aavalid_op_sub:
	case aavalid_op_mul:
	case aavalid_op_div:
	case aavalid_op_max:
	case aavalid_op_min:
        case aavalid_op_cancel2:
          fprintf(stderr, "x_fix = "); ia_print(stderr, d->x_fix); fprintf(stderr,"\n");
          fprintf(stderr, "y_fix = "); ia_print(stderr, d->y_fix); fprintf(stderr,"\n");
	  break;
	case aavalid_op_affine:
	case aavalid_op_neg:
	case aavalid_op_inv:
	case aavalid_op_sqr:
	case aavalid_op_sqrt:
	case aavalid_op_abs:
          fprintf(stderr, "x_fix = "); ia_print(stderr, d->x_fix); fprintf(stderr,"\n");
	  break;
	case aavalid_op_cancel1:
	  fprintf(stderr, "y_fix = "); ia_print(stderr, d->y_fix); fprintf(stderr, "\n");
	  break;
	default:
	  fatalerror("aavalid_throw_and_fix_eps: bad op code");
      }
    fprintf(stderr, "z_cmp = "); ia_print(stderr, d->z_cmp); fprintf(stderr, "\n");
    fprintf(stderr, "z_fix = "); ia_print(stderr, d->z_fix); fprintf(stderr, "\n");
  }

int aavalid_is_zero(AAP x)
  {
    return(aa_is_zero(x)); 
  }
    
Float aavalid_throw_float(void)
  { 
    switch (rand()&3)
      { case 0: return(-One);
        case 1: return(Zero);
        case 2: return(One);
        case 3: return(flt_random_mag(0, 5) * (flt_random() - flt_random()));
        default: return(Zero);
      }
  }
  
Float aavalid_throw_eps(void)
  {
    double n;

    ROUND_NEAR;
    n = (rand()&65535) / 65536.0;
    if ((n -= 0.25) < 0.0) 
      return(-1);
    else if ((n -= 0.25) < 0.0)
      return(0);
    else if ((n -= 0.25) < 0.0)
      return(1);
    else
      { Float num = (Float)((rand() & 8388607) / 8388608.0);
	if (n < 0.125) num = -num;
	return(num);
      }
  }





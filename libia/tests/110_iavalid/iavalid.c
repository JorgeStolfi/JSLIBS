/* Validation of IA ops */
/* Last edited on 2024-12-21 11:23:58 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <flt.h>
#include <ia.h>
#include <affirm.h>
#include <jsrandom.h>

#define DEBUG 0

typedef enum {
    iavalid_op_const =         0,
    iavalid_op_int_const =     1,
    iavalid_op_add =           2,
    iavalid_op_sub =           3,
    iavalid_op_neg =           4,
    iavalid_op_scale =         5,
    iavalid_op_shift =         6,
    iavalid_op_affine =        7,
    iavalid_op_affine_2 =      8,
    iavalid_op_inv =           9,
    iavalid_op_mul =          10,
    iavalid_op_div =          11,
    iavalid_op_sqr =          12,
    iavalid_op_sqrt =         13,
    iavalid_op_exp =          14,
    iavalid_op_abs =          15,
    iavalid_op_max =          16,
    iavalid_op_min =          17,
    iavalid_op_interp =       18 
  } iavalid_op_t;

typedef struct IAData {  /* Data for one IA operands/result tuple */
    iavalid_op_t  op;  /* Operation code. */
    char *op_name;   /* Operation's name. */
    /* Arguments for {ia_const}: */
    Float fltx;      /* Value argument for {ia_const}. */
    Float flte;      /* Error argument for {ia_const}. */
    /* Arguments for {ia_int_const}: */
    int intx;        /* Argument for {ia_int_const}. */
    /* Arguments for {ia_scale,ia_shift,ia_affine,ia_affine_2}: */
    Float alpha;     /* Multiplier for {x}. */
    Float beta;      /* Multiplier for {y}. */
    Float zeta;      /* Denominator. */
    Float gamma;     /* Constant term. */
    Float delta;     /* Extra error term. */
    /* Arguments for {ia_interp}: */
    Float x0;        /* Abscissa of first interval {y0} (here {x}). */
    Float x1;        /* Abscissa of second interval {y1} (here {y}). */
    Float xev;       /* Interpolation abscissa {x}. */
    /* Main arguments and results: */
    Interval x;           /* First argument of {op}. */
    Interval y;           /* Second argument of {op} (if binary). */
    Interval z;           /* Result of {ia_op(x, y)}. */
    /* Validation intervals: */
    Interval x_fix;  /* Sample value in {x}. */
    Interval y_fix;  /* Sample value in {y}. */
    Interval z_cmp;  /* Result of {ia_op(x_fix, y_fix)}. */
  } IAData;

/*** INTERNAL PROTOTYPES ***/

int main(int argc, char *argv[]);

void iavalid_test_op (
    iavalid_op_t op, /* The operation to test */
    char *op_name,   /* The operation's name */
    int nforms,      /* Number of affine evaluations */
    int npoints      /* Number of sample points */
  );
  /* Tests an IA operation {nforms} times.  For each test, generates
     a set {x, y, ...} of IA operands, evaluates {op} on them, and then
     performs {npoints} pointwise checks on the resulting IA form {z}.
     
     Each pointwise check consists of selecting a random {Float} value
     in each operand {x, y, ..}; converting those {Float} values to IA
     representation, which results in almost point-like intervals
     {x_fix, y_fix, ...}; applying {op} to those intervals, in IA
     arithmetic, which results in an interval {z_cmp}; and finally
     checking that {z_cmp} and {z} have nonempty intersection.
     
     Note that this is not a strong test, since it is likely to
     succeed even if the IA operation is completely --- but
     consistently --- wrong. */

void iavalid_throw_and_eval_ia_tuple(IAData *d);
  /* Generates an adequate IA operand tuple {d->x}, {d->y} 
    for {d->op}, and performs {ia_op} on it, storing the result
    in {d->z}. */
    
void iavalid_sample_intervals(IAData *d);
  /* Selects random {Float} values in {d->x} and {d->y},
    obtaining {d->x_fix} and {d->y_fix}.  Then
    performs {d->op} with IA arithmetic on {d->x_fix} and
    {d->y_fix}, obtaining {d->z_cmp}. */
    
void iavalid_print_data(IAData *d);
  /* Prints {*d} to {stderr}. */
    
Interval iavalid_sample(Interval x);
  /* Generates a random pointlike interval contained in {x}. */
    
Float iavalid_throw_float(void);
  /* Generates a random float that is 1, 0, -1, or anything, with
    probability 0.25 each. */

int iavalid_is_zero(Interval x);
  /* TRUE if the IA {x} is exactly zero. */

/*** MAIN PROGRAM ***/

int main(int argc, char *argv[])
  {
    flt_init();
    ia_init();
  
    iavalid_test_op(iavalid_op_const,        "const",     100,    1);
    iavalid_test_op(iavalid_op_int_const,    "int_const", 100,    1);
    iavalid_test_op(iavalid_op_add,          "add",       200, 5000);
    iavalid_test_op(iavalid_op_sub,          "sub",       200, 5000);
    iavalid_test_op(iavalid_op_neg,          "neg",       200, 5000);
    iavalid_test_op(iavalid_op_shift,        "shift",     200, 5000);
    iavalid_test_op(iavalid_op_scale,        "scale",     200, 5000);
    iavalid_test_op(iavalid_op_affine,       "affine",    200, 5000);
    iavalid_test_op(iavalid_op_affine_2,     "affine_2",  200, 5000);
    iavalid_test_op(iavalid_op_inv,          "inv",       200, 5000);
    iavalid_test_op(iavalid_op_mul,          "mul",       200, 5000);
    iavalid_test_op(iavalid_op_div,          "div",       200, 5000);
    iavalid_test_op(iavalid_op_sqr,          "sqr",       200, 5000);
    iavalid_test_op(iavalid_op_sqrt,         "sqrt",      200, 5000);
    iavalid_test_op(iavalid_op_sqrt,         "exp",       200, 5000);
    iavalid_test_op(iavalid_op_abs,          "abs",       200, 5000);
    iavalid_test_op(iavalid_op_max,          "max",       200, 5000);
    iavalid_test_op(iavalid_op_min,          "min",       200, 5000);
    iavalid_test_op(iavalid_op_interp,       "interp",    200, 5000);
    
    return(0);

  }
		     
void iavalid_test_op (
    iavalid_op_t op,   /* The operation to test */
    char *op_name,   /* The operation's name */
    int nforms,      /* Number of affine evalutaions */
    int npoints      /* Number of sample points */
  )
  {
    IAData d;
    int iform, ipoint;
    int bad, debug;

    fprintf(stderr, "\n");
    fprintf(stderr, "--------------------------------------------------------\n");
    fprintf(stderr, "testing %s\n", op_name);
    fprintf(stderr, "\n");
    
    srand(314159);
    
    d.op = op;
    d.op_name = op_name;
    
    for (iform=0; iform < nforms; iform++)
      {
        iavalid_throw_and_eval_ia_tuple(&d);
        for (ipoint=0; ipoint < npoints; ipoint++)
          { 
            iavalid_sample_intervals(&d);
            bad = ((d.z_cmp.hi < d.z.lo) || (d.z_cmp.lo > d.z.hi));
            debug = ((ipoint == 0) & (iform < 3));
            if (bad || debug) 
              { fprintf(stderr, "iform = %d  ipoint = %d\n", iform, ipoint);
                iavalid_print_data(&d);
                affirm(!bad, "z and z_cmp are disjoint");
                fprintf(stderr, "\n");
              }
          }
      }

    fprintf(stderr, "\n");
    fprintf(stderr, "--------------------------------------------------------\n");
  }
  
void iavalid_throw_and_eval_ia_tuple(IAData *d)
  {
    switch (d->op)
      {
	case iavalid_op_const:
	  d->fltx = iavalid_throw_float();
          d->flte = fabsf(iavalid_throw_float());
          d->z = ia_const(d->fltx, d->flte);
	  break;
	case iavalid_op_int_const:
	  d->intx = rand()-rand();
          d->z = ia_int_const(d->intx);
	  break;
	case iavalid_op_scale:
	  d->alpha = iavalid_throw_float();
          do { d->zeta = iavalid_throw_float(); } while (d->zeta == Zero);
          d->x = ia_throw();
          d->z = ia_scale(d->x, d->alpha, d->zeta);
	  break;
	case iavalid_op_shift:
	  d->gamma = iavalid_throw_float();
          d->x = ia_throw();
          d->z = ia_shift(d->x, d->gamma);
	  break;
	case iavalid_op_affine:
	  d->alpha = iavalid_throw_float();
          do { d->zeta = iavalid_throw_float(); } while (d->zeta == Zero);
          d->gamma = iavalid_throw_float();
          d->delta = iavalid_throw_float();
          d->x = ia_throw();
          d->z = ia_affine(d->x, d->alpha, d->zeta, d->gamma, d->delta);
	  break;
	case iavalid_op_affine_2:
	  d->alpha = iavalid_throw_float();
          d->beta = iavalid_throw_float();
          do { d->zeta = iavalid_throw_float(); } while (d->zeta == Zero);
          d->gamma = iavalid_throw_float();
          d->delta = iavalid_throw_float();
          d->x = ia_throw();
          d->y = ia_throw();
          d->z = ia_affine_2(d->x, d->alpha, d->y, d->beta, d->zeta, d->gamma, d->delta);
	  break;
	case iavalid_op_add:
	  d->x = ia_throw();
          d->y = ia_throw();
          d->z = ia_add(d->x, d->y);
	  break;
	case iavalid_op_sub:
	  d->x = ia_throw();
          d->y = ia_throw();
          d->z = ia_sub(d->x, d->y);
	  break;
	case iavalid_op_neg:
	  d->x = ia_throw();
          d->z = ia_neg(d->x);
	  break;
	case iavalid_op_inv:
          /* The {x} interval must contain some non-zero value: */
	  do { d->x = ia_throw(); } while (iavalid_is_zero(d->x));
          d->z = ia_inv(d->x);
	  break;
	case iavalid_op_mul:
	  d->x = ia_throw();
          d->y = ia_throw();
          d->z = ia_mul(d->x, d->y);
	  break;
	case iavalid_op_div:
	  d->x = ia_throw();
          /* The {y} interval must contain some non-zero value: */
	  do { d->y = ia_throw(); } while (iavalid_is_zero(d->y));
          d->z = ia_div(d->x, d->y);
	  break;
	case iavalid_op_sqr:
	  d->x = ia_throw();
          d->z = ia_sqr(d->x);
	  break;
	case iavalid_op_sqrt:
	  /* The {x} interval must contain some non-negative value: */
	  do { d->x = ia_throw(); } while (d->x.hi < Zero);
          d->z = ia_sqrt(d->x);
	  break;
	case iavalid_op_exp:
	  d->x = ia_throw();
          d->z = ia_exp(d->x);
	  break;
	case iavalid_op_abs:
	  d->x = ia_throw();
          d->z = ia_abs(d->x);
	  break;
	case iavalid_op_max:
	  d->x = ia_throw();
          d->y = ia_throw();
          d->z = ia_max(d->x, d->y);
	  break;
	case iavalid_op_min:
	  d->x = ia_throw();
          d->y = ia_throw();
          d->z = ia_min(d->x, d->y);
	  break;
	case iavalid_op_interp:
	  d->x0 = iavalid_throw_float();
          d->x1 = iavalid_throw_float();
          d->xev = iavalid_throw_float(); /* Alias {x} of {ia_interp}. */
          d->x = ia_throw(); /* Alias {y0} of {ia_interp}. */
          d->y = ia_throw(); /* Alias {x0} of {ia_interp}. */
          d->z = ia_interp(d->x0, d->x, d->x1, d->y, d->xev);
	  break;
	default:
	  fatalerror("iavalid_throw_and_eval_ia_tuple: bad op code");
      }
  }

void iavalid_sample_intervals(IAData *d)
  {
    switch (d->op)
      {
	case iavalid_op_const:
	  { Interval xv = (Interval){ d->fltx, d->fltx }; 
            Interval ev = (Interval){ -d->flte, +d->flte };
            Interval dv = iavalid_sample(ev);
            d->z_cmp = ia_add(xv, dv);
          }
	  break;
	case iavalid_op_int_const:
          { Float xi = (Float)(d->intx);
            d->z_cmp = (Interval){ xi, xi };
          }
	  break;
	case iavalid_op_shift:
	  { d->x_fix = iavalid_sample(d->x);
            d->z_cmp = ia_shift(d->x_fix, d->gamma);
          }
	  break;
	case iavalid_op_scale:
	  { d->x_fix = iavalid_sample(d->x);
            d->z_cmp = ia_scale(d->x_fix, d->alpha, d->zeta);
	  }
          break;
	case iavalid_op_affine:
	  { d->x_fix = iavalid_sample(d->x);
            d->z_cmp = ia_affine(d->x_fix, d->alpha, d->zeta, d->gamma, d->delta);
	  }
          break;
	case iavalid_op_affine_2:
	  { d->x_fix = iavalid_sample(d->x);
            d->y_fix = iavalid_sample(d->y);
            d->z_cmp = ia_affine_2(d->x_fix, d->alpha, d->y_fix, d->beta, d->zeta, d->gamma, d->delta);
	  }
          break;
	case iavalid_op_add:
	  { d->x_fix = iavalid_sample(d->x);
            d->y_fix = iavalid_sample(d->y);
            d->z_cmp = ia_add(d->x_fix, d->y_fix);
	  }
          break;
	case iavalid_op_sub:
	  { d->x_fix = iavalid_sample(d->x);
            d->y_fix = iavalid_sample(d->y);
            d->z_cmp = ia_sub(d->x_fix, d->y_fix);
	  }
          break;
	case iavalid_op_neg:
	  { d->x_fix = iavalid_sample(d->x);
            d->z_cmp = ia_neg(d->x_fix);
	  }
          break;
	case iavalid_op_inv:
	  { do { d->x_fix = iavalid_sample(d->x); } while (iavalid_is_zero(d->x_fix)); 
            d->z_cmp = ia_inv(d->x_fix);
	  }
          break;
	case iavalid_op_mul:
	  { d->x_fix = iavalid_sample(d->x);
            d->y_fix = iavalid_sample(d->y);
            d->z_cmp = ia_mul(d->x_fix, d->y_fix);
	  }
          break;
	case iavalid_op_div:
	  { d->x_fix = iavalid_sample(d->x);
            do { d->y_fix = iavalid_sample(d->y); } while (iavalid_is_zero(d->y_fix)); 
            d->z_cmp = ia_div(d->x_fix, d->y_fix);
	  }
          break;
	case iavalid_op_sqr:
	  { d->x_fix = iavalid_sample(d->x);
            d->z_cmp = ia_sqr(d->x_fix);
	  }
          break;
	case iavalid_op_sqrt:
	  { /* Remove the negative part of the range before sampling: */
            Interval tx = d->x;
            if (tx.lo < 0) { tx.lo = 0; }
            do { d->x_fix = iavalid_sample(tx); } while (d->x_fix.hi < Zero); 
            d->z_cmp = ia_sqrt(d->x_fix);
	  }
          break;
	case iavalid_op_exp:
	  { d->x_fix = iavalid_sample(d->x);
            d->z_cmp = ia_exp(d->x_fix);
	  }
          break;
	case iavalid_op_abs:
	  { d->x_fix = iavalid_sample(d->x);
            d->z_cmp = ia_abs(d->x_fix);
	  }
          break;
	case iavalid_op_max:
	  { d->x_fix = iavalid_sample(d->x);
            d->y_fix = iavalid_sample(d->y);
            d->z_cmp = ia_max(d->x_fix, d->y_fix);
	  }
          break;
	case iavalid_op_min:
	  { d->x_fix = iavalid_sample(d->x);
            d->y_fix = iavalid_sample(d->y);
            d->z_cmp = ia_min(d->x_fix, d->y_fix);
	  }
          break;
	case iavalid_op_interp:
	  { d->x_fix = iavalid_sample(d->x);
            d->y_fix = iavalid_sample(d->y);
            d->z_cmp = ia_interp(d->x0, d->x_fix, d->x1, d->y_fix, d->xev);
	  }
          break;
	default:
	  fatalerror("iavalid_sample_intervals: bad op code");
      }
  }

void iavalid_print_data(IAData *d)
  {
    switch (d->op)
      {
	case iavalid_op_const:
	  fprintf(stderr, "fltx = "); flt_print(stderr, d->fltx); fprintf(stderr, "\n");
	  fprintf(stderr, "flte = "); flt_print(stderr, d->flte); fprintf(stderr, "\n");
	  break;
	case iavalid_op_int_const:
          fprintf(stderr, "intx = %d\n", d->intx);
	  break;
	case iavalid_op_shift:
	  fprintf(stderr, "gamma = "); flt_print(stderr, d->gamma); fprintf(stderr, "\n");
	  fprintf(stderr, "x = "); ia_print(stderr, d->x); fprintf(stderr, "\n");
	  break;
	case iavalid_op_scale:
	  fprintf(stderr, "alpha = "); flt_print(stderr, d->alpha); fprintf(stderr, "\n");
	  fprintf(stderr, "zeta  = "); flt_print(stderr, d->zeta);  fprintf(stderr, "\n");
	  fprintf(stderr, "x = "); ia_print(stderr, d->x); fprintf(stderr, "\n");
	  break;
	case iavalid_op_affine:
	  fprintf(stderr, "alpha = "); flt_print(stderr, d->alpha); fprintf(stderr, "\n");
	  fprintf(stderr, "zeta  = "); flt_print(stderr, d->zeta);  fprintf(stderr, "\n");
	  fprintf(stderr, "gamma = "); flt_print(stderr, d->gamma); fprintf(stderr, "\n");
	  fprintf(stderr, "delta = "); flt_print(stderr, d->delta); fprintf(stderr, "\n");
          fprintf(stderr, "x = "); ia_print(stderr, d->x); fprintf(stderr, "\n");
	  break;
	case iavalid_op_affine_2:
	  fprintf(stderr, "alpha = "); flt_print(stderr, d->alpha); fprintf(stderr, "\n");
	  fprintf(stderr, "beta  = "); flt_print(stderr, d->beta);  fprintf(stderr, "\n");
	  fprintf(stderr, "zeta  = "); flt_print(stderr, d->zeta);  fprintf(stderr, "\n");
	  fprintf(stderr, "gamma = "); flt_print(stderr, d->gamma); fprintf(stderr, "\n");
	  fprintf(stderr, "delta = "); flt_print(stderr, d->delta); fprintf(stderr, "\n");
          fprintf(stderr, "x = "); ia_print(stderr,d->x); fprintf(stderr,"\n");
          fprintf(stderr, "y = "); ia_print(stderr,d->y); fprintf(stderr,"\n");
	  break;
	case iavalid_op_add:
	case iavalid_op_sub:
	case iavalid_op_mul:
	case iavalid_op_div:
	case iavalid_op_max:
	case iavalid_op_min:
          fprintf(stderr, "x = "); ia_print(stderr, d->x); fprintf(stderr, "\n");
          fprintf(stderr, "y = "); ia_print(stderr, d->y); fprintf(stderr, "\n");
	  break;
	case iavalid_op_neg:
	case iavalid_op_inv:
	case iavalid_op_sqr:
	case iavalid_op_sqrt:
	case iavalid_op_exp:
	case iavalid_op_abs:
          fprintf(stderr, "x = "); ia_print(stderr,d->x); fprintf(stderr,"\n");
	  break;
	case iavalid_op_interp:
	  fprintf(stderr, "x0 = "); flt_print(stderr, d->x0);  fprintf(stderr, "\n");
	  fprintf(stderr, "y0 = "); ia_print(stderr, d->x);  fprintf(stderr, "\n");
	  fprintf(stderr, "x1 = "); flt_print(stderr, d->x1);  fprintf(stderr, "\n");
	  fprintf(stderr, "y1 = "); ia_print(stderr, d->y); fprintf(stderr, "\n");
	  fprintf(stderr, "x =  "); flt_print(stderr, d->xev); fprintf(stderr, "\n");
	  break;
	default:
	  fatalerror("iavalid_print_intervals: bad op code");
      }
    fprintf(stderr,"z = "); ia_print(stderr,d->z); fprintf(stderr,"\n");
    fprintf(stderr, "\n");
    switch (d->op)
      {
	case iavalid_op_const:
	case iavalid_op_int_const:
	  break;
	case iavalid_op_affine_2:
	case iavalid_op_add:
	case iavalid_op_sub:
	case iavalid_op_mul:
	case iavalid_op_div:
	case iavalid_op_max:
	case iavalid_op_min:
        case iavalid_op_interp:
          fprintf(stderr, "x_fix = "); ia_print(stderr, d->x_fix); fprintf(stderr,"\n");
          fprintf(stderr, "y_fix = "); ia_print(stderr, d->y_fix); fprintf(stderr,"\n");
	  break;
	case iavalid_op_affine:
        case iavalid_op_shift:
        case iavalid_op_scale:
	case iavalid_op_neg:
	case iavalid_op_inv:
	case iavalid_op_sqr:
	case iavalid_op_sqrt:
	case iavalid_op_abs:
          fprintf(stderr, "x_fix = "); ia_print(stderr, d->x_fix); fprintf(stderr,"\n");
	  break;
	default:
	  fatalerror("iavalid_print_intervals: bad op code");
      }
    fprintf(stderr, "z_cmp = "); ia_print(stderr, d->z_cmp); fprintf(stderr, "\n");
  }

int iavalid_is_zero(Interval x)
  {
    return (x.lo == Zero) && (x.hi == Zero); 
  }
    
Float iavalid_throw_float(void)
  { 
    switch (rand()&3)
      { case 0: return(-One);
        case 1: return(Zero);
        case 2: return(One);
        case 3: return(flt_random_mag(0, 5) * (flt_random() - flt_random()));
        default: return(Zero);
      }
  }
  
Interval iavalid_sample(Interval x)
  {
    demand(x.lo <= x.hi, "bad interval");
    if (x.lo == x.hi) { return x; }
    Float rad = ia_rad(x);
    Float mid = ia_mid(x);
    Float del;
    switch (rand()&3)
      { case 0: del = -One;
        case 1: del = Zero;
        case 2: del = +One;
        case 3: del = flt_random() - flt_random();
        default: del = Zero;
      }
    Float smp = mid + del*rad;
    if (smp < x.lo) { smp = x.lo; }
    if (smp > x.hi) { smp = x.hi; }
    return (Interval){ smp, smp };
  }





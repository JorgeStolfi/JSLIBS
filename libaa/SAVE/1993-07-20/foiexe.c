/* FOI exercise program */

#include "foi.h"
#include "foifloat.h"
#include "foimisc.h"
#include "interval.h"
#include "iomisc.h"
#include <math.h>
#include <stdio.h>

#define NTERMS 6

typedef struct {Float eps[NTERMS]; } Assignment; 
  /* Assignment of values to the noise vars [0..NTERMS-1], in [-1__+1] */

/*** INTERNAL PROTOTYPES ***/

int main(void);
void foix_error(char *msg);
void foix_test_zero(int effort);
void foix_test_bin(
    char *op, 
    int effort, 
    FOIP ff (FOIP x, FOIP y), 
    Interval fv (Interval, Interval y)
  );

void foix_test_un(
    char *op, 
    int effort, 
    FOIP ff (FOIP x), 
    Interval fv (Interval)
  );

void foix_assign_eps(Assignment *asgn);
Interval foix_fix_foi(FOIP x, Assignment *asgn);

/*** IMPLEMENTATIONS ***/

int main(void)
  {
    char op[20];
    int effort;

    srandom(2201);

    foi_init();
    scanf("%[a-z] d%, op, &effort);
    if (strcmp(op, "zero")  == 0) foix_test_zero (effort);
    if (strcmp(op, "add")   == 0) foix_test_bin  (op, effort, foi_add,  iv_add);
    if (strcmp(op, "sub")   == 0) foix_test_bin  (op, effort, foi_sub,  iv_sub);
    if (strcmp(op, "neg")   == 0) foix_test_un   (op, effort, foi_neg,  iv_neg);
    if (strcmp(op, "mul")   == 0) foix_test_bin  (op, effort, foi_mul,  iv_mul);
    if (strcmp(op, "sqr")   == 0) foix_test_un   (op, effort, foi_sqr,  iv_sqr);
    if (strcmp(op, "inv")   == 0) foix_test_un   (op, effort, foi_inv,  iv_inv);
    if (strcmp(op, "sqrt")  == 0) foix_test_un   (op, effort, foi_sqrt, iv_sqrt);

    fclose(stdout);
    return (0);
  }

void foix_error(char *msg)
  {
    printf("\n*** error -- ");
    printf(msg);
    printf(" ***\n");
    fflush(stdout);
    *((int*)-1) = *((int*)-1); /* bomb out */
    exit(1);
  }

void foix_test_zero(int effort)
  {
    Interval xv;
    FOIP x, y, z;

    printf("\n=== Testing foi_zero ===\n\n");

    x = foi_zero();
    printf("foi_zero() = "); foi_print(stdout, x); putchar('\n');
    xv = foi_implicit_range();
    if((xv.lo != Zero)|(xv.hi != Zero)) foix_error("foi_zero: bad range");
  }
  
void foix_test_bin(
    char *op, 
    int effort, 
    FOIP ff (FOIP x, FOIP y),
    Interval fv (Interval x, Interval y)
  )
  {
    Interval xv, yv, zv, zr;
    Assignment asgn;
    FOIP x, y, z;
    int i;
    frame: MemP;

    printf("\n=== Testing foi_%s ===\n\n", op);

    for (i=0; i<effort; i++)
      {
	frame = foi_top();
        x = foix_throw(NTERMS);
	y = foix_throw(NTERMS);
	z = ff (x, y);
        zr = foi_range(z);
        for (j=0; j<5000; j++)
          {
            foix_assign_eps(&asgn);
            xv = foix_fix_foi(x, &asgn);
            yv = foix_fix_foi(y, &asgn);
            zv = fv(xv, yv);
            if ((zr.lo > zv.hi)|(zr.hi < zv.lo)) foix_error("inconsistency");
          }
        foi_flush(frame);
      }
  }
  
void foix_test_un(
    char *op, 
    int effort, 
    FOIP ff (FOIP x),
    Interval fv (Interval x)
  )
  {
    Interval xv, zv, zr;
    Assignment asgn;
    FOIP x, z;
    int i;
    frame: MemP;

    printf("\n=== Testing foi_%s ===\n\n", op);

    for (i=0; i<effort; i++)
      {
	frame = foi_top();
        x = foix_throw(NTERMS);
	z = ff (x);
        zr = foi_range(z);
        for (j=0; j<5000; j++)
          {
            foix_assign_eps(&asgn);
            xv = foix_fix_foi(x, &asgn);
            zv = fv(xv);
            if ((zr.lo > zv.hi)|(zr.hi < zv.lo)) foix_error("inconsistency");
          }
        foi_flush(frame);
      }
  }
  
void foix_assign_eps(Assignment *asgn)
  {
    for (i=0; i<NTERMS; i++) 
      asgn->eps[i] = foix_throw_eps();
  }
  
Float foix_throw_eps(void)
  {
    switch (random()&3)
      {
        case 0: return (-1.0);
        case 1: return (+1.0);
        case 2: return (0.0);
        case 3: return (flt_random());
      }
  }

Interval foix_fix_foi(FOIP x, Assignment *asgn)
  {
    Interval res = iv_const(x->center, Zero);
    TermP xp = ((TermP)(x+1));
    TermCount xn = x->nterms;
    while (xn>0)
      {
        if (xp->id >= NTERMS) foix_error("foix_fix_foi: bad eps id");
        res = iv_add(
          res,
          iv_mul(
            iv_const(xp->coef, 0.0), 
            iv_const(asgn->eps[xp->id], 0.0)
          )
        );
        xp++;
        xn--;
      }
  }






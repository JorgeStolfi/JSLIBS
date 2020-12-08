/* Last edited on 2012-07-21 12:28:31 by stolfi */
/* Random simple tests of the AA library */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <flt.h>
#include <ia.h>
#include <aa.h>

/*** INTERNAL PROTOTYPES ***/

int main(int argc, char *argv[]);

void aa_test_inv(void);
void aa_test_div(void);
void aa_test_return_n(void);

void aa_test_print(char *name, AAP x);

/*** MAIN PROGRAM ***/

int main(int argc, char *argv[])
  {
    flt_init();
    ia_init();
    aa_init();
  
    aa_test_return_n();
    aa_test_inv();
    aa_test_div();
    fprintf(stderr,"----------------------------------------------------------------------\n");
    return (0);
  }

void aa_test_inv(void)
  {
    fprintf(stderr,"--- aa_test_inv ------------------------------------------------------\n");
    MemP frame = aa_top();
    Interval xi = (Interval){3,4};
    AAP xa = aa_from_interval(xi); 
    AAP za;
    za = aa_inv(xa);
    aa_test_print("x", xa);
    aa_test_print("z", za);
    aa_flush(frame);
    fprintf(stderr,"----------------------------------------------------------------------\n");
  }

void aa_test_div(void)
  {
    fprintf(stderr,"--- aa_test_div ------------------------------------------------------\n");
    MemP frame = aa_top();
    Interval xi = (Interval){1,1};
    Interval yi = (Interval){3,4};
    AAP xa = aa_from_interval(xi);
    AAP ya = aa_from_interval(yi);
    AAP za;
    za = aa_div(xa,ya);
    aa_test_print("x", xa);
    aa_test_print("y", ya);
    aa_test_print("z", za);
    aa_flush(frame);
  }

void aa_test_return_n(void)
  {
    fprintf(stderr,"--- aa_test_return_n -------------------------------------------------\n");
    MemP frame = aa_top();
    fprintf(stderr,"frame = %012lu\n", (uintptr_t)frame);
    /* Put some garbage on the stack: */
    Interval xi = (Interval){1,1};
    Interval yi = (Interval){3,4};
    AAP xa = aa_from_interval(xi);
    AAP ya = aa_from_interval(yi);
    AAP za = aa_mul(xa, ya);
    aa_test_print("x", xa);
    aa_test_print("y", ya);
    aa_test_print("z", za);
    fprintf(stderr,"\n");
    /* Remember the current top-of-stack: */
    MemP subframe = aa_top();
    fprintf(stderr,"subframe = %012lu\n", (uintptr_t)frame);
    /* Do some computations: */
    AAP t[6];
    t[0] = aa_add(xa,ya);
    t[1] = aa_sub(za,xa);
    t[2] = aa_affine(xa, Two, One, Zero, Zero);
    t[3] = aa_add(t[2],t[0]);
    t[4] = aa_affine(ya, Four, One, Zero, Zero);
    t[5] = aa_mul(t[1],t[4]);
    /* Return some intersting forms, in the wrong order: */
    int n = 2;
    AAP r[n];
    r[0] = t[4];
    r[1] = t[2];
    aa_test_print("r[0]", r[0]);
    aa_test_print("r[1]", r[1]);
    fprintf(stderr,"\n");
    fprintf(stderr,"calling aa_return_n ...\n");
    aa_return_n(subframe,n,r);
    fprintf(stderr,"\n");
    /* See whether the forms have survived: */
    aa_test_print("x", xa);
    aa_test_print("y", ya);
    aa_test_print("z", za);
    fprintf(stderr,"\n");
    aa_test_print("r[0]", r[0]);
    aa_test_print("r[1]", r[1]);
    aa_flush(frame);
    fprintf(stderr,"----------------------------------------------------------------------\n");
  }

void aa_test_print(char *name, AAP x)
  {
    uintptr_t xp = (uintptr_t)x;
    fprintf(stderr, "%-5s (%012lu) = ", name, xp);
    aa_print(stderr, x); 
    fprintf(stderr, "\n");
  }

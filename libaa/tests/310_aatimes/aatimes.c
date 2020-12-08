/* Timing tests for libaa (the affine arithmetic library) */
/* Last edited on 2016-12-26 22:10:37 by stolfilocal */

#include <stdio.h>
#include <stdint.h>

#include <flt.h>
#include <timefunc.h>
#include <aa.h>

/* #include <sys/ieeefp.h> */

int ieee_flags(char *a, char *b, char *c, char **out);

int main (void);

int64_t do_empty(int64_t ntimes);
int64_t do_aa_add(int64_t ntimes);
int64_t do_aa_mul(int64_t ntimes);

#define NTIMES 200000
#define NEPS 10

int main (void)
  {
    aa_init();
    
    time_func("empty loop",     do_empty,  do_empty,     NTIMES);
    time_func("aa_add",         do_aa_add, do_empty,     NTIMES);
    time_func("aa_mul",         do_aa_mul, do_empty,     NTIMES);
    return(0);
  }
    
int64_t do_empty(int64_t ntimes)
  { int64_t i;
    MemP frame = aa_top();
    AAP x1 = aa_throw(NEPS);
    AAP x2 = aa_throw(NEPS);
    AAP x3 = x2;
    AAP x4 = x1;
    MemP locframe = aa_top();
    AAP z;
    for (i=0; i<ntimes; i++) 
      { z = x1; 
        x1 = x2; 
        x2 = x3;
        x3 = x4; 
        x4 = z; 
        z = x1; 
        x1 = x2; 
        x2 = x3;
        x3 = x4; 
        x4 = z; 
        aa_flush(locframe);
      }
    aa_flush(frame);
    return (10*ntimes);
  }
 
int64_t do_aa_add(int64_t ntimes)
  { int64_t i;
    MemP frame = aa_top();
    AAP x = aa_throw(NEPS);
    AAP y = aa_throw(NEPS);
    MemP locframe = aa_top();
    for (i=0; i<ntimes; i++)
      { (void)aa_add(x, y);
        (void)aa_add(x, y);
        (void)aa_add(x, y);
        (void)aa_add(x, y);
        (void)aa_add(x, y);
        (void)aa_add(x, y);
        (void)aa_add(x, y);
        (void)aa_add(x, y);
        (void)aa_add(x, y);
        (void)aa_add(x, y);
        aa_flush(locframe);
      }
    aa_flush(frame);
    return (10*ntimes);
  }
    
int64_t do_aa_mul(int64_t ntimes)
  { int64_t i;
    MemP frame = aa_top();
    AAP x = aa_throw(NEPS);
    AAP y = aa_throw(NEPS);
    MemP locframe = aa_top();
    for (i=0; i<ntimes; i++)
      { (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        (void)aa_mul(x, y);
        aa_flush(locframe);
      }
    aa_flush(frame);
    return (10*ntimes);
  }
    


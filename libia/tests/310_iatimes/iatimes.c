/* iatimes.c - timing of IA operations */
/* Last edited on 2024-12-21 11:24:05 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <ia.h>
#include <flt.h>
#include <jsrandom.h>

#include <timefunc.h>

/* int ieee_flags(char *a, char *b, char *c, char **out); */

int main (void);

int64_t do_empty(int64_t ntimes);
int64_t do_ia_add(int64_t ntimes);
int64_t do_ia_mul(int64_t ntimes);

#define NTIMES 10000000

int main (void)
  {
    ia_init();
    
    time_func("empty loop",     do_empty,  do_empty,     NTIMES);
    time_func("ia_add",         do_ia_add, do_empty,     NTIMES);
    time_func("ia_mul",         do_ia_mul, do_empty,     NTIMES);
    return(0);
  }
    
int64_t do_empty(int64_t ntimes)
  { int64_t i;
    Interval x1 = (Interval){flt_random()-0.5f, flt_random()+0.5f};
    Interval x2 = (Interval){flt_random()-0.5f, flt_random()+0.5f};
    Interval x3 = x1;
    Interval x4 = x2;
    Interval z;
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
      }
    return (10*ntimes);
  }
 
int64_t do_ia_add(int64_t ntimes)
  { int64_t i;
    Interval x = (Interval){flt_random()-0.5f, flt_random()+0.5f};
    Interval y = (Interval){flt_random()-0.5f, flt_random()+0.5f};
    Interval z;
    for (i=0; i<ntimes; i++)
      { z = ia_add(x, y);
        z = ia_add(x, z);
        z = ia_add(x, z);
        z = ia_add(x, z);
        z = ia_add(x, z);
        z = ia_add(x, z);
        z = ia_add(x, z);
        z = ia_add(x, z);
        z = ia_add(x, z);
        (void)ia_add(x, z);
      }
    return (10*ntimes);
  }
    
int64_t do_ia_mul(int64_t ntimes)
  { int64_t i;
    Interval x = (Interval){flt_random()-0.5f, flt_random()+0.5f};
    Interval y = (Interval){flt_random()-0.5f, flt_random()+0.5f};
    Interval z;
    for (i=0; i<ntimes; i++)
      { z = ia_mul(x, y);
        z = ia_mul(x, z);
        z = ia_mul(x, z);
        z = ia_mul(x, z);
        z = ia_mul(x, z);
        z = ia_mul(x, z);
        z = ia_mul(x, z);
        z = ia_mul(x, z);
        z = ia_mul(x, z);
        (void)ia_mul(x, z);
      }
    return (10*ntimes);
  }
    


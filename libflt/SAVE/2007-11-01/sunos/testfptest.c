/* Quick test of rounding mode setting (fptest.s prototype) - Sun SPARC4/Solaris5 */
/* Last edited on 2007-01-04 00:14:47 by stolfi */

#include <ioprotos.h>
#include <stdio.h>
#include <ieeefp.h>
#include <sys/ieeefp.h>
#include <flt.h>
#include "fptest.h"

int main(int argc, char **argv);

int main(int argc, char **argv)
{ float x, y, z, one, three;
  int fsr;
  one = 1.0;
  three = 3.0;
  
  /* Testing the Sun routines: */

  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  fpsetround(FP_RP);         fprintf(stderr, "fpsetround(FP_RP)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  fpsetround(FP_RM);         fprintf(stderr, "fpsetround(FP_RM)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  fpsetround(FP_RZ);         fprintf(stderr, "fpsetround(FP_RZ)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  fpsetround(FP_RN);         fprintf(stderr, "fpsetround(FP_RN)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 

  /* Testing the "fptest.s" routines: */

  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  set_rounding(fp_positive); fprintf(stderr, "set_rounding(fp_positive)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  set_rounding(fp_negative); fprintf(stderr, "set_rounding(fp_negative)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  set_rounding(fp_tozero);   fprintf(stderr, "set_rounding(fp_tozero)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  set_rounding(fp_nearest);  fprintf(stderr, "set_rounding(fp_nearest)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 

  /* Testing the "flt.h" routines: */

  flt_init();
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  ROUND_UP;                  fprintf(stderr, "ROUND_UP\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  ROUND_DOWN;                fprintf(stderr, "ROUND_DOWN\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  ROUND_ZERO;                fprintf(stderr, "ROUND_ZERO\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  ROUND_NEAR;                fprintf(stderr, "ROUND_NEAR\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 

  set_rounding(fp_positive);
  x = one/three;
  set_rounding(fp_negative);
  y = one/three;
  set_rounding(fp_nearest);
  z = x - y;
  
  fprintf(stderr, "up   = %18.8e\n", x); 
  fprintf(stderr, "down = %18.8e\n", y); 
  fprintf(stderr, "diff = %18.8e\n", z);
  fclose(stderr);
  return (0);
}

 

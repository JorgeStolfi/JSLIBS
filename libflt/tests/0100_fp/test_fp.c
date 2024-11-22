/* Quick test of rounding mode setting (fptest.s prototype) - Intel/Linux */
/* Last edited on 2024-11-20 07:50:59 by stolfi */

#include <stdio.h>
#include <fenv.h>
#include <fpu_control.h>
#include <flt.h>
#include <test_fp_tools.h>

int main(int argc, char **argv);

int main(int argc, char **argv)
{ float x, y, z, one, three;
  int fsr;
  one = 1.0;
  three = 3.0;
  
  /* Testing the Linux rounding routines: */

  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  fesetround(FE_UPWARD);     fprintf(stderr, "fesetround(FE_UPWARD)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  fesetround(FE_DOWNWARD);   fprintf(stderr, "fesetround(FE_DOWNWARD)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  fesetround(FE_TOWARDZERO); fprintf(stderr, "fesetround(FE_TOWARDZERO)\n"); 
  fsr = get_fsr();           fprintf(stderr, "fsr   = %8x\n", fsr); 
  fesetround(FE_TONEAREST);  fprintf(stderr, "fesetround(FE_TONEAREST)\n"); 
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

  set_rounding(FE_UPWARD);
  x = one/three;
  set_rounding(FE_DOWNWARD);
  y = one/three;
  set_rounding(FE_TONEAREST);
  z = x - y;
  
  fprintf(stderr, "up   = %18.8e\n", x); 
  fprintf(stderr, "down = %18.8e\n", y); 
  fprintf(stderr, "diff = %18.8e\n", z);
  fclose(stderr);
  return (0);
}

 

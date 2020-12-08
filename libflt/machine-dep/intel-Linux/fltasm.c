/* Setting rounding-mode direction */
/* Last edited on 2005-02-14 23:33:13 by stolfi */

#include <flt.h>
#include <fenv.h>
#include <fpu_control.h>

void flt_round_down (void)
{ 
  fesetround(FE_DOWNWARD);
}
  
void flt_round_up   (void)
{ 
  fesetround(FE_UPWARD);
}
  
void flt_round_near (void)
{ 
  fesetround(FE_TONEAREST);
}
  
void flt_round_zero (void)
{ 
  fesetround(FE_TOWARDZERO);
}

int flt_get_fsr(void)
{
  fpu_control_t cw;
  _FPU_GETCW(cw);
  return cw;
}

void flt_single_prec (void)
{
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw = (cw & ~_FPU_EXTENDED & ~_FPU_DOUBLE) | _FPU_SINGLE;
  _FPU_SETCW(cw);
}

void flt_double_prec (void)
{
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw = (cw & ~_FPU_EXTENDED & ~_FPU_SINGLE) | _FPU_DOUBLE;
  _FPU_SETCW(cw);
}

void flt_extended_prec (void)
{
  fpu_control_t cw;
  _FPU_GETCW(cw);
  cw = (cw & ~_FPU_DOUBLE & ~_FPU_SINGLE) | _FPU_EXTENDED;
  _FPU_SETCW(cw);
}

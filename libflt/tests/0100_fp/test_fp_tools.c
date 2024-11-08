/* Setting rounding-mode direction */
/* Last edited on 2024-11-08 12:31:32 by stolfi */

#include <fenv.h>
#include <fpu_control.h>
#include <test_fp_tools.h>

void set_rounding (int dir)
{ 
  fesetround(dir);
}

int get_fsr(void)
{
  fpu_control_t cw;
  _FPU_GETCW(cw);
  return cw;
}

/* basic floating-point definitions - Sun SPARC4/Solaris 5 */
/* Last edited on 2024-11-08 12:31:19 by stolfi */

#ifndef test_fp_tools_H
#define test_fp_tools_H

#include <fpu_control.h>

void set_rounding (int dir);
int get_fsr (void); 

#endif

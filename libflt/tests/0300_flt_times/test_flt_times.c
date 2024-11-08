/* Timing of floating point ops (from flt.h) */
/* Last edited on 2016-12-26 18:30:26 by stolfilocal */

#include <stdio.h>
#include <stdint.h>

#include <timefunc.h>

#include <flt.h>

// #if (defined(SunOS4))
// #  include <sys/ieeefp.h>
// int ieee_flags(char *a, char *b, char *c, char **out);
// #endif
// 
// #if (defined(SunOS5))
// #  include <ieeefp.h>
// #endif

#include <fenv.h>
#include <fpu_control.h>

#define APPSECS 5 
  /* Approximate time for each test in seconds. */

// #if (defined(OSF1V4))
// /* #  include <ieeefp.h> */
// #endif

int main (int argc, char **argv);

int64_t do_empty(int64_t nloop);

int64_t do_conversion(int64_t nloop);
int64_t do_single_add(int64_t nloop);
int64_t do_single_mul(int64_t nloop);
int64_t do_single_div(int64_t nloop);
int64_t do_double_add(int64_t nloop);
int64_t do_double_mul(int64_t nloop);
int64_t do_double_div(int64_t nloop);
int64_t do_FMIN_FMAX(int64_t nloop);
int64_t do_my_set_round(int64_t nloop);
int64_t do_sys_set_round(int64_t nloop);
int64_t do_flt_mul(int64_t nloop);
int64_t do_flt_mix(int64_t nloop);

int main (int argc, char **argv)
  {
    flt_init();
    
    { /* Tests flt_mul overflow: */
      Float x, y, z, err;
      x = (Float)(MaxFloat/2.0);
      y = (Float)(MaxFloat/4.0);
      err = 0;
      flt_mul(x, y, &z, &err);
      ROUND_NEAR;
      fprintf(stderr, " x =    %18.8f\n", x);
      fprintf(stderr, " y =    %18.8f\n", y);
      fprintf(stderr, " z =    %18.8f\n", z);
      fprintf(stderr, " err =  %18.8f\n", err);
      fprintf(stderr, " +Inf = %18.8f\n", PlusInfinity);
      fprintf(stderr, " -Inf = %18.8f\n", MinusInfinity);
      fprintf(stderr, " MaxF = %18.8f\n", MaxFloat);
    }

    time_func("empty loop",     do_empty,         do_empty,   75000000*APPSECS);
    time_func("conversion",     do_conversion,    do_empty,   50000000*APPSECS);
    time_func("single add",     do_single_add,    do_empty,   10000000*APPSECS);
    time_func("single mul",     do_single_mul,    do_empty,   10000000*APPSECS);
    time_func("single div",     do_single_div,    do_empty,    5000000*APPSECS);
    time_func("double add",     do_double_add,    do_empty,   10000000*APPSECS);
    time_func("double mul",     do_double_mul,    do_empty,   10000000*APPSECS);
    time_func("double div",     do_double_div,    do_empty,    5000000*APPSECS);
    time_func("FMIN/FMAX",      do_FMIN_FMAX,     do_empty,   50000000*APPSECS);
    time_func("ROUND_UP/DOWN",  do_my_set_round,  do_empty,    1500000*APPSECS);
    time_func("sys set_round",  do_sys_set_round, do_empty,    1500000*APPSECS);
    time_func("flt_mul",        do_flt_mul,       do_empty,     300000*APPSECS);
    time_func("flt_mix",        do_flt_mix,       do_empty,     500000*APPSECS);
    return(0);
  }
    
int64_t do_empty(int64_t nloop)
  { int64_t i;
    Float x1, x2, x3, x4, x5, x6, x7, x8, x9, t;
    x1 = (Float)0.1;
    x2 = (Float)0.2;
    x3 = (Float)0.3;
    x4 = (Float)0.4;
    x5 = (Float)0.5;
    x6 = (Float)0.6;
    x7 = (Float)0.7;
    x8 = (Float)0.8;
    x9 = (Float)0.9;
    for (i=0; i<nloop; i++) {
      t = x1;
      x1 = x2;
      x2 = x3;
      x3 = x4;
      x4 = x5;
      x5 = x6;
      x6 = x7;
      x7 = x8;
      x8 = x9;
      x9 = t;
    }
    return (10*nloop);
  }
 
int64_t do_conversion(int64_t nloop)
  { int64_t i;
    float x;
    double xx = 0.7;
    for (i=0; i<nloop; i++)
      { x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
        x = (float) xx; xx = (double) x;
      }
    return (10*nloop);
  }
    
int64_t do_single_add(int64_t nloop)
  { int64_t i;
    Float x, y;
    x = (Float)0.1;
    y = (Float)0.1;
    for (i=0; i<nloop; i++)
      { x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
      }
    return (10*nloop);
  }
    
int64_t do_single_mul(int64_t nloop)
  { int64_t i;
    Float x, y1, y2;
    x = (Float)0.1;
    y1 = (Float)1.001;
    y2 = (Float)(1.0/y1);
    for (i=0; i<nloop; i++)
      { x = x * y1;
        x = x * y2;
        x = x * y1;
        x = x * y2;
        x = x * y1;
        x = x * y2;
        x = x * y1;
        x = x * y2;
        x = x * y1;
        x = x * y2;
      }
    return (10*nloop);
  }
    
int64_t do_single_div(int64_t nloop)
  { int64_t i;
    Float x, y1, y2;
    x = (Float)0.1;
    y1 = (Float)1.001;
    y2 = (Float)(1.0/y1);
    for (i=0; i<nloop; i++)
      { x = x / y1;
        x = x / y2;
        x = x / y1;
        x = x / y2;
        x = x / y1;
        x = x / y2;
        x = x / y1;
        x = x / y2;
        x = x / y1;
        x = x / y2;
      }
    return (10*nloop);
  }
    
int64_t do_double_add(int64_t nloop)
  { int64_t i;
    double x, y;
    x = 0.1;
    y = 0.1;
    for (i=0; i<nloop; i++)
      { x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
        x = x + y;
      }
    return (10*nloop);
  }
    
int64_t do_double_mul(int64_t nloop)
  { int64_t i;
    double x, y1, y2;
    x = 0.1;
    y1 = 1.001;
    y2 = 1.0/y1;
    for (i=0; i<nloop; i++)
      { x = x * y1;
        x = x * y2;
        x = x * y1;
        x = x * y2;
        x = x * y1;
        x = x * y2;
        x = x * y1;
        x = x * y2;
        x = x * y1;
        x = x * y2;
      }
    return (10*nloop);
  }
    
int64_t do_double_div(int64_t nloop)
  { int64_t i;
    double x, y1, y2;
    x = 0.1;
    y1 = 1.001;
    y2 = 1.0/y1;
    for (i=0; i<nloop; i++)
      { x = x / y1;
        x = x / y2;
        x = x / y1;
        x = x / y2;
        x = x / y1;
        x = x / y2;
        x = x / y1;
        x = x / y2;
        x = x / y1;
        x = x / y2;
      }
    return (10*nloop);
  }

int64_t do_FMIN_FMAX(int64_t nloop)
  { int64_t i;
    Float x, y;
    x = (Float)0.7;
    y = (Float)0.3;
    for (i=0; i<nloop; i++)
      { x = FMIN(x, y);
        x = FMAX(x, y);
        x = FMIN(x, y);
        x = FMAX(x, y);
        x = FMIN(x, y);
        x = FMAX(x, y);
        x = FMIN(x, y);
        x = FMAX(x, y);
        x = FMIN(x, y);
        x = FMAX(x, y);
      }
    return (10*nloop);
  }
    
int64_t do_my_set_round(int64_t nloop)
  { int64_t i;
    for (i=0; i<nloop; i++)
      { ROUND_UP;
        ROUND_DOWN;
        ROUND_UP;
        ROUND_DOWN;
        ROUND_UP;
        ROUND_DOWN;
        ROUND_UP;
        ROUND_DOWN;
        ROUND_UP;
        ROUND_DOWN;
      }
    return (10*nloop);
  }

// #if (defined(SunOS4))
// int64_t do_sys_set_round(int64_t nloop)
//   { int64_t i;
//     char *out_dir;
//     for (i=0; i<nloop; i++)
//       { 
//         ieee_flags("set", "direction", "positive", &out_dir);
//         ieee_flags("set", "direction", "negative", &out_dir);
// 
//         ieee_flags("set", "direction", "positive", &out_dir);
//         ieee_flags("set", "direction", "negative", &out_dir);
// 
//         ieee_flags("set", "direction", "positive", &out_dir);
//         ieee_flags("set", "direction", "negative", &out_dir);
// 
//         ieee_flags("set", "direction", "positive", &out_dir);
//         ieee_flags("set", "direction", "negative", &out_dir);
// 
//         ieee_flags("set", "direction", "positive", &out_dir);
//         ieee_flags("set", "direction", "negative", &out_dir);
//       }
//     return (10*nloop);
//   }
// #endif
// 
// #if (defined(SunOS5))
// int64_t do_sys_set_round(int64_t nloop)
//   { int64_t i;
//     fp_rnd out_dir;
//     for (i=0; i<nloop; i++)
//       { 
//         out_dir = fpsetround(FP_RP);
//         out_dir = fpsetround(FP_RN);
//         out_dir = fpsetround(FP_RZ);
//         out_dir = fpsetround(FP_RM);
//         out_dir = fpsetround(FP_RZ);
// 
//         out_dir = fpsetround(FP_RP);
//         out_dir = fpsetround(FP_RN);
//         out_dir = fpsetround(FP_RZ);
//         out_dir = fpsetround(FP_RM);
//         out_dir = fpsetround(FP_RZ);
//       }
//     return (10*nloop);
//   }
// #endif

int64_t do_sys_set_round(int64_t nloop)
  { int64_t i;
    int out_dir[10];
    for (i=0; i<nloop; i++)
      { 
        out_dir[0] = fesetround(FE_UPWARD);
        out_dir[1] = fesetround(FE_DOWNWARD);
        out_dir[2] = fesetround(FE_UPWARD);
        out_dir[3] = fesetround(FE_DOWNWARD);
        out_dir[4] = fesetround(FE_UPWARD);
        out_dir[5] = fesetround(FE_DOWNWARD);
        out_dir[6] = fesetround(FE_UPWARD);
        out_dir[7] = fesetround(FE_DOWNWARD);
        out_dir[8] = fesetround(FE_UPWARD);
        out_dir[9] = fesetround(FE_DOWNWARD);
      }
    out_dir[0] = out_dir[9]; /* To pacify the compiler. */
    return (10*nloop);
  }

int64_t do_flt_mul(int64_t nloop)
  { int64_t i;
    Float x, y1, y2, err;
    x = (Float)0.1;
    y1 = (Float)1.001;
    y2 = (Float)(1.0/y1);
    err = (Float)0.0;
    for (i=0; i<nloop; i++)
      { flt_mul(x, y1, &x, &err);
        flt_mul(x, y2, &x, &err);
        flt_mul(x, y1, &x, &err);
        flt_mul(x, y2, &x, &err);
        flt_mul(x, y1, &x, &err);
        flt_mul(x, y2, &x, &err);
        flt_mul(x, y1, &x, &err);
        flt_mul(x, y2, &x, &err);
        flt_mul(x, y1, &x, &err);
        flt_mul(x, y2, &x, &err);
      }
    return (10*nloop);
  }

int64_t do_flt_mix(int64_t nloop)
  { int64_t i;
    Float x, y, a1, a2, b1, b2, z1, z2, err;
    x = (Float)0.1;
    y = (Float)0.9;
    
    a1 = (Float)+0.625;
    b1 = (Float)+0.375;
    z1 = (Float)+3.000;
    
    a2 = (Float)+3.000;
    b2 = (Float)-0.375;
    z2 = (Float)+0.625;

    err = (Float)0.0;
    for (i=0; i<nloop; i++)
      { flt_mix(x, a1, y, b1, z1, &x, &err);
        flt_mix(x, a2, y, b2, z2, &x, &err);
        flt_mix(x, a1, y, b1, z1, &x, &err);
        flt_mix(x, a2, y, b2, z2, &x, &err);
        flt_mix(x, a1, y, b1, z1, &x, &err);
        flt_mix(x, a2, y, b2, z2, &x, &err);
        flt_mix(x, a1, y, b1, z1, &x, &err);
        flt_mix(x, a2, y, b2, z2, &x, &err);
        flt_mix(x, a1, y, b1, z1, &x, &err);
        flt_mix(x, a2, y, b2, z2, &x, &err);
      }
    return (10*nloop);
  }


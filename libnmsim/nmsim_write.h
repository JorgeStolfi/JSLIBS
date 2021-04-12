#ifndef nmsim_write_H
#define nmsim_write_H
 
/* Value writing procedures. */
/* Last edited on 2020-12-16 22:47:41 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#define nmsim_write_VIJ_PREC 0.0006
  /* Max rounding error allowed for {V,I,J,I_avg,I_dev} values. */

#define nmsim_write_age_PREC 0.06
  /* Max rounding error allowed for mean {age} values. */

#define nmsim_write_MH_PREC 0.0000006
  /* Max rounding error allowed for {M,H,M_avg,H_avg} values. */

#define nmsim_write_tau_PREC 0.0000006
  /* Max rounding error allowed for {tau} values. */

#define nmsim_write_rho_PREC 0.000006
  /* Max rounding error allowed for {rho,eta} values. */

#define nmsim_write_KW_PREC 0.000006
  /* Max rounding error allowed for {K,W,Wavg,Wdev} values. */

void nmsim_write_double_value(FILE *wr, double v, double prec, bool_t sgn, bool_t fudge_0, bool_t fudge_1);
  /* Writes {v} using format "%.{N}f", where {N} is such that the rounding error (0.5 ulp) 
    is no more than {prec}. Eliminates trailing fraction
    zeros. If all fraction digits are blanked out, also blanks the
    decimal point.
    
    If the absolute value of {v}, after formatting and blanking,
    prints as "0", no sign is written. Otherwise, writes a "-" sign
    only if {v} is negative, and a "+" only if {v} is positive and
    {sgn} is true.
  
    If {fudge_0} is true, writes the exact value 0 as "0", and 
    prints any other values that would print as "0" as
    the nearest nonzero value, preserving their signs.  If {fudge_0}
    is false, may write non-zero values as "0".
    
    If {fudge_1} is true, writes the exact value 1 and as "1" (or
    "+1"), and prints any other values that would print as "1" (or
    "+1") as the nearest value different from 1. If {fudge_0} is
    false, may write values that are different from 1 as "1" (or
    "+1"). Ditto for -1. */

void nmsim_write_int64_param(FILE *wr, char *pref, char *name, int64_t v, char *fmt);
  /* Writes a line "{pref}{name} = {v}" using the given format (of the 'd' kind). */

void nmsim_write_double_param
  ( FILE *wr, 
    char *pref, 
    char *name,
    double v, 
    double prec,
    bool_t sgn, 
    bool_t fudge_0, 
    bool_t fudge_1
  );
  /* Writes a line "{pref}{name} = {v}" using 
    {nmsim_write_double_value(wr,v,prec,sgn,fudge_0,fudge_1)}.
    Removes leading blanks and trailing fractional zeros. */

#endif

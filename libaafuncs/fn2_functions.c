/* See {fn2_functions.h}. */
/* Last edited on 2007-01-04 03:03:36 by stolfi */

#include <fn2_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <affirm.h>
#include <bool.h>

#include <fn2_f0.h>
#include <fn2_f1.h>
#include <fn2_f2.h>
#include <fn2_f3.h>
#include <fn2_f4.h>
#include <fn2_f5.h>
#include <fn2_f6.h>
#include <fn2_f7.h>
#include <fn2_f8.h>
#include <fn2_fadd.h>
#include <fn2_fdiv.h>
#include <fn2_fiamany.h>
#include <fn2_fmax.h>
#include <fn2_fmax2.h>
#include <fn2_fmul.h>
#include <fn2_fmul2.h>

#include <string.h>

fn2_data_t fn2_from_tag(char *tag)
  {
    if (strcmp(tag, "f0")      == 0) { return fn2_f0_get_data(); }
    if (strcmp(tag, "f1")      == 0) { return fn2_f1_get_data(); }
    if (strcmp(tag, "f2")      == 0) { return fn2_f2_get_data(); }
    if (strcmp(tag, "f3")      == 0) { return fn2_f3_get_data(); }
    if (strcmp(tag, "f4")      == 0) { return fn2_f4_get_data(); }
    if (strcmp(tag, "f5")      == 0) { return fn2_f5_get_data(); }
    if (strcmp(tag, "f6")      == 0) { return fn2_f6_get_data(); }
    if (strcmp(tag, "f7")      == 0) { return fn2_f7_get_data(); }
    if (strcmp(tag, "f8")      == 0) { return fn2_f8_get_data(); }
    if (strcmp(tag, "fadd")    == 0) { return fn2_fadd_get_data(); }
    if (strcmp(tag, "fdiv")    == 0) { return fn2_fdiv_get_data(); }
    if (strcmp(tag, "fiamany") == 0) { return fn2_fiamany_get_data(); }
    if (strcmp(tag, "fmax")    == 0) { return fn2_fmax_get_data(); }
    if (strcmp(tag, "fmax2")   == 0) { return fn2_fmax2_get_data(); }
    if (strcmp(tag, "fmul")    == 0) { return fn2_fmul_get_data(); }
    if (strcmp(tag, "fmul2")   == 0) { return fn2_fmul2_get_data(); }
    demand(FALSE, "bad function tag");
    return fn2_f0_get_data(); /* To pacify the compiler. */
  }

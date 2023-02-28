/* See bbgoal.h */
/* Last edited on 2023-02-20 06:50:40 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <string.h>

#include <bool.h>
#include <affirm.h>

#include <fbb_f3_ia.h>
#include <fbb_f2_ia.h>
#include <fbb_f1_ia.h>

#include <bbgoal.h>

bbgoal_data_t bbgoal_from_tag(char *tag)
  {
    if (strcmp(tag, "f1_ia") == 0) { return fbb_f1_ia_get_data(); }
    if (strcmp(tag, "f2_ia") == 0) { return fbb_f2_ia_get_data(); }
    if (strcmp(tag, "f3_ia") == 0) { return fbb_f3_ia_get_data(); }
    demand(FALSE, "bad function tag");
    return fbb_f1_ia_get_data(); /* To pacify the compiler. */
  }

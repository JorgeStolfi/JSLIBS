/* See {image_output_coords.h}. */
/* Last edited on 2023-08-28 14:40:57 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <argparser.h>

#include <image_output_coords.h>

#define SMA  (1.0e-100)
  /* A very small value, but still far from underflow. */
  
#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

void imgc_parse_output_unit(argparser_t *pp, double *oUnit)
  { 
    if (argparser_keyword_present(pp, "-oUnit"))
      { (*oUnit) = argparser_get_next_double(pp, SMA, BIG); }
    else
      { (*oUnit) = 1.0; }
  }

void imgc_parse_output_center_org(argparser_t *pp, bool_t *oCenter, r2_t *oOrg)
  {
    if (argparser_keyword_present(pp, "-oCenter"))
      { (*oCenter) = TRUE; }
    else if (argparser_keyword_present(pp, "-oOrg"))
      { (*oCenter) = FALSE;
        oOrg->c[0] = argparser_get_next_double(pp, -BIG, +BIG);
        oOrg->c[1] = argparser_get_next_double(pp, -BIG, +BIG);
      }
  }

void imgc_parse_output_size(argparser_t *pp, double *oCols, double *oRows)
  {
    if (argparser_keyword_present(pp, "-oSize"))
      { (*oCols) = argparser_get_next_double(pp, SMA, BIG);
        (*oRows) = argparser_get_next_double(pp, SMA, BIG);
      }
  }

/* See {image_input_coords.h}. */
/* Last edited on 2023-08-28 14:40:32 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <argparser.h>
#include <image_input_coords.h>

#define SMA  (1.0e-100)
  /* A very small value, but still far from underflow. */
  
#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

void imgc_parse_input_unit(argparser_t *pp, double *iUnit)
  { 
    if (argparser_keyword_present(pp, "-iUnit"))
      { (*iUnit) = argparser_get_next_double(pp, SMA, BIG); }
    else
      { (*iUnit) = 1.0; }
  }

void imgc_parse_input_center_org(argparser_t *pp, bool_t *iCenter, r2_t *iOrg)
  {
    if (argparser_keyword_present(pp, "-iCenter"))
      { (*iCenter) = TRUE; }
    else if (argparser_keyword_present(pp, "-iOrg"))
      { (*iCenter) = FALSE;
        iOrg->c[0] = argparser_get_next_double(pp, -BIG, +BIG);
        iOrg->c[1] = argparser_get_next_double(pp, -BIG, +BIG);
      }
  }

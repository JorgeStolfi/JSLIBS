/* See {image_coords.h}. */
/* Last edited on 2023-08-27 17:33:04 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <r3x3.h>
#include <hr2.h>
#include <argparser.h>
#include <image_coords.h>

#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

#define imgc_unit_MAX (1024.0*1024.0*1024.0)
#define imgc_unit_MIN (1.0/imgc_unit_MAX)
  /* Paranoia bounds for "-iUnit", "-oUnit". */

hr2_pmap_t imgc_coord_sys_map
  ( bool_t xRev, 
    bool_t yRev, 
    double unit,
    bool_t center, 
    r2_t *org, 
    int cols, 
    int rows
  )
  {
    /* Start with the identity matrix {A}, get its weight {w}: */
    r3x3_t A; r3x3_ident(&A);
    assert(A.c[0][0] == 1.0);
      
    for (int32_t ax = 1; ax <= 2; ax++)
      { bool_t rev = (ax == 1 ? xRev : yRev);
        double oc = org->c[ax-1];            /* Position of origin in pixel units. */
        double sz = (ax == 1 ? cols : rows); /* Image size in pixels. */
        double d = (rev ? -1/unit : +1/unit);
        double t; /* Translation term. */
        A.c[ax][ax] = d;
        if (center)
          { t = oc - d*0.5*sz; }
        else
          { t = (rev ? oc - d*sz : oc); }
        A.c[0][ax] = t;
      }
   
    /* Build the map: */
    hr2_pmap_t M;
    M.dir = A; r3x3_inv(&A, &(M.inv));
    return M;
  }

void imgc_parse_x_axis(argparser_t *pp, bool_t *xLeft)
  {
    if 
      ( (argparser_keyword_present(pp, "-xAxis")) ||
        (argparser_keyword_present(pp, "-hAxis"))
      )
      { if (argparser_keyword_present_next(pp, "left"))
          { (*xLeft) = TRUE; }
        else if (argparser_keyword_present_next(pp, "right"))
          { (*xLeft) = FALSE; }
        else
          { argparser_error(pp, "invalid horizontal axis direction"); }
      }
  }
  
void imgc_parse_y_axis(argparser_t *pp, bool_t *yDown)
  {
    if  
      ( (argparser_keyword_present(pp, "-yAxis"))||
        (argparser_keyword_present(pp, "-vAxis"))
      )
      { if (argparser_keyword_present_next(pp, "up"))
          { (*yDown) = FALSE; }
        else if (argparser_keyword_present_next(pp, "down"))
          { (*yDown) = TRUE; }
        else
          { argparser_error(pp, "invalid vertical axis direction"); }
      }
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

void imgc_parse_output_size(argparser_t *pp, int *oCols, int *oRows, int max_size)
  {
    if (argparser_keyword_present(pp, "-oSize"))
      { (*oCols) = (int)argparser_get_next_int(pp, 1, max_size);
        (*oRows) = (int)argparser_get_next_int(pp, 1, max_size);
      }
  }

void imgc_parse_input_unit(argparser_t *pp, double *iUnit)
  { 
    if (argparser_keyword_present(pp, "-iUnit"))
      { (*iUnit) = argparser_get_next_double(pp, imgc_unit_MIN, imgc_unit_MAX); }
    else
      { (*iUnit) = 1.0; }
  }

void imgc_parse_output_unit(argparser_t *pp, double *oUnit)
  { 
    if (argparser_keyword_present(pp, "-oUnit"))
      { (*oUnit) = argparser_get_next_double(pp, imgc_unit_MIN, imgc_unit_MAX); }
    else
      { (*oUnit) = 1.0; }
  }

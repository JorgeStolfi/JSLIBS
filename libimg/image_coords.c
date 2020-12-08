/* See {image_coords.h}. */
/* Last edited on 2013-10-21 02:44:11 by stolfilocal */

#include <bool.h>
#include <r2.h>
#include <r3x3.h>
#include <hr2.h>
#include <argparser.h>
#include <image_coords.h>

#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

hr2_pmap_t imgc_coord_sys_map
  ( bool_t xRev, 
    bool_t yRev, 
    bool_t center, 
    r2_t *org, 
    int cols, 
    int rows
  )
  {
    /* Start with the identity matrix {A}, get its weight {w}: */
    r3x3_t A; r3x3_ident(&A);
    double w = A.c[0][0];
    
    /* Apply the axis-reversal and default origin changes: */
    if (xRev) { A.c[1][1] = -w; A.c[0][1] = w*cols; }
    if (yRev) { A.c[2][2] = -w; A.c[0][2] = w*rows; }
    
    /* Apply the origin shift: */
    if (center) 
      { A.c[0][1] -= 0.5*w*cols; A.c[0][2] -= 0.5*w*rows; }
    else
      { A.c[0][1] -= w*org->c[0]; A.c[0][2] -= w*org->c[1]; }
      
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

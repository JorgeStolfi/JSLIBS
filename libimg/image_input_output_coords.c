/* See {image_input_output_coords.h}. */
/* Last edited on 2023-09-25 09:22:20 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <r3x3.h>
#include <hr2.h>
#include <argparser.h>
#include <image_input_output_coords.h>
#include <image_input_coords.h>
#include <image_output_coords.h>

#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

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

void imgc_compute_output_size_in_pixels
  ( int32_t iCols_pix,
    int32_t iRows_pix,
    double iUnit,
    double oCols_usr,
    double oRows_usr,
    double oUnit,
    int32_t *oCols_pix_P,
    int32_t *oRows_pix_P,
    int32_t max_pix
  )
  { 
    /* Input image size in user units: */
    double iCols_usr = iCols_pix/iUnit;
    double iRows_usr = iRows_pix/iUnit;
    
    /* Provide default size of output image (in user units) if not specified: */
    if (oCols_usr < 0) { oCols_usr = iCols_usr; }
    if (oRows_usr < 0) { oRows_usr = iRows_usr; }
    
    /* Compute output image size in pixels, trying to ignore roundoff errors: */
    int32_t oCols_pix = (int32_t)fmax(1.0, ceil(oCols_usr*oUnit - 0.0001));
    int32_t oRows_pix = (int32_t)fmax(1.0, ceil(oRows_usr*oUnit - 0.0001));
    
    /* Range check: */
    demand(oCols_pix <= max_pix, "too many output pixel columns"); 
    demand(oRows_pix <= max_pix, "too many output pixel rows"); 
    
    (*oCols_pix_P) = oCols_pix;
    (*oRows_pix_P) = oRows_pix;
  }

void imgc_print_matrix(FILE *wr, char *name1, char *name2, char *name3, char *tag, r3x3_t *M)
  { 
    fprintf(stderr, "  %s to %s coordinates (%s%s) =\n", name1, name2, name3, tag);
    r3x3_gen_print
      ( wr, M, 
        "%10.4f",
        "", "", "",
        "    [ ", " ", " ]\n" 
      );
    fprintf(stderr, "\n");
  }

void imgc_print_pmap(FILE *wr, char *name1, char *name2, char *name3, hr2_pmap_t *M)
  { fprintf(wr, "\n");
    imgc_print_matrix(wr, name1, name2, name3, "", &(M->dir)); 
    imgc_print_matrix(wr, name2, name1, name3, "^-1", &(M->inv));
  }
    
void imgc_parse_input_output_coords_args
  ( argparser_t *pp,
    bool_t *xLeft,
    bool_t *yUp,
    double *iUnit,
    bool_t *iCenter,
    r2_t *iOrg,
    double *oUnit,
    bool_t *oCenter,
    r2_t *oOrg,
    double *oCols,
    double *oRows
  )
  {
    imgc_parse_x_axis(pp, xLeft);
    imgc_parse_y_axis(pp, yUp);
    imgc_parse_input_unit(pp, iUnit);
    imgc_parse_input_center_org(pp, iCenter, iOrg);
    imgc_parse_output_unit(pp, oUnit);
    imgc_parse_output_center_org(pp, oCenter, oOrg);
    imgc_parse_output_size(pp, oCols, oRows);
  }

void imgc_parse_x_axis(argparser_t *pp, bool_t *xLeft)
  {
    if
      ( argparser_keyword_present(pp, "-xAxis") ||
        argparser_keyword_present(pp, "-hAxis")
      )
      { if
          ( argparser_keyword_present_next(pp, "left") ||
            argparser_keyword_present_next(pp, "l")
          )
          { (*xLeft) = TRUE; }
        else if 
          (  argparser_keyword_present_next(pp, "right")||
            argparser_keyword_present_next(pp, "r")
          )
          { (*xLeft) = FALSE; }
        else
          { argparser_error(pp, "invalid horizontal axis direction"); }
      }
  }
  
void imgc_parse_y_axis(argparser_t *pp, bool_t *yUp)
  {
    if  
      ( argparser_keyword_present(pp, "-yAxis") ||
        argparser_keyword_present(pp, "-vAxis")
      )
      { if 
          ( argparser_keyword_present_next(pp, "up") ||
            argparser_keyword_present_next(pp, "u")
          )
          { (*yUp) = TRUE; }
        else if 
          ( argparser_keyword_present_next(pp, "down")||
            argparser_keyword_present_next(pp, "d")
          )
          { (*yUp) = FALSE; }
        else
          { argparser_error(pp, "invalid vertical axis direction"); }
      }
  }


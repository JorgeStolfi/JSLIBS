/* See pst_lamp.h */
/* Last edited on 2025-01-03 12:38:13 by stolfi */

#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h> 
#include <r3.h> 
#include <frgb.h>
#include <frgb_ops.h>
#include <affirm.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_lamp.h>
#include <pst_argparser.h>

vec_typeimpl(pst_lamp_vec_t, pst_lamp_vec, pst_lamp_t *);

pst_lamp_t *pst_lamp_new(r3_t *dir, frgb_t *pwr, double crad)
  { pst_lamp_t *src = talloc(1, pst_lamp_t);
    
    if (dir == NULL) 
      { src->dir = (r3_t){{ NAN, NAN, NAN }}; }
    else
      { src->dir = *(dir); 
        double dm = r3_dir(&(src->dir), &(src->dir ));
        if ((! isfinite(dm)) || (dm == 0)) { src->dir = (r3_t){{ NAN, NAN, NAN }}; }
      }
    
    if (pwr == NULL)
      { src->pwr = frgb_NoColor; }
    else
      { src->pwr = (*pwr);
        if (! (isfinite(src->pwr.c[0]) && isfinite(src->pwr.c[1]) && isfinite(src->pwr.c[2])))
          { src->pwr = frgb_NoColor; }
      }
    src->crad = crad;
    demand(isnan(crad) || ((crad >= -1) && (crad <= +1)), "invalid {crad}");
    return src;
  }

double pst_lamp_geom_factor(r3_t *nrm, r3_t *dir, double crad)
  { double coef; /* Geometric illumination coefficient. */
    if (crad <= -1.0) 
      { /* Lamp spans the whole sphere (ambient term). */ 
        coef = 1.0;
      }
    else 
      { double clum = r3_dot(nrm, dir); /* Cosine of illumination angle. */
        if (crad == 1.0)
          { /* point source: */
            return (clum > 0 ? clum : 0.0);
          }
        else if (crad == 0.0)
          { /* Lamp spans a whole hemisphere (wall term). */
            return 0.25*(1 + clum)*(1 + clum);
          }
        else if (crad > 0.0)
          { /* Lamp spans at most one hemisphere: */
            double srad = sqrt(fmax(0, 1-crad*crad)); /* Sine of ang radius of source. */
            if (clum < -srad)
              { /* Point is in full shadow: */
                coef = 0.0;
              }
            else if (clum > srad)
              { /* Point sees the whole lamp: */
                coef = clum; /* By definition, {coef=1} if {clum=1}. */
              }
            else
              { /* Penumbra - a C1 transition between {clum=-srad} and {clum=+srad}: */
                double s = clum/srad;
                coef = 0.25*srad*(1 + s)*(1 + s);
              }
          }
        else
          { /* Lamp spans more than one hemisphere: */
            double srad = sqrt(fmax(0, 1-crad*crad)); /* Sine of ang radius of source. */
            if (clum > srad)
              { /* Point sees a full hemisphere of lamp; by definition {coef=1}. */
                coef = 1.0; 
              }
            else if (clum < -srad)
              { /* Point sees the full "hole" of the lamp: */
                coef = 1 + clum*(1 + crad);
              }
            else
              { 
                double s = 1.0 - clum/srad;
                coef = 1.0 - (1 + crad)*0.25*srad*(1 + s)*(1 + s);
              }
          }
      }
    return coef;
  }
  
#define pst_bogus_spec_MESS \
  " is not applicable or was already specified for this lamp"

pst_lamp_t *pst_lamp_spec_parse(argparser_t *pp, uint32_t *N_P)
  { 
    pst_lamp_t *src = pst_lamp_new(NULL, NULL, NAN);
    src->dir = (r3_t){{ NAN, NAN, NAN }};
    src->pwr = (frgb_t){{ NAN, NAN, NAN }};

    /* Parse the lamp type: */
    bool_t dir_given = FALSE;
    bool_t crad_given = FALSE;
    bool_t pwr_given = FALSE;

    uint32_t N = 1;
    if (argparser_keyword_present_next(pp, "lamp"))
      { /* Nothing to do. */
        src->crad = NAN;
      }
    else if ((N_P != NULL) && argparser_keyword_present_next(pp, "array"))
      { N = (uint32_t)argparser_get_next_int(pp, 1, UINT32_MAX); }
    else if (argparser_keyword_present_next(pp, "wall"))
      { src->crad = 0.0; }
    else if (argparser_keyword_present_next(pp, "ambient"))
      { src->crad = -1.0;
        src->dir = (r3_t){{ 0.0, 0.0, 0.0 }};
      }
    else
      { free(src); src = NULL; N = 0; }
    
    if (src != NULL)
      { 
        while (TRUE)
          { /* One more keyword: */
            if (argparser_keyword_present_next(pp, "angles"))
              { if (dir_given) { argparser_error(pp, "lamp's direction" pst_bogus_spec_MESS); }
                /* Get azimuth and elevation in degrees: */
                double azim = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
                double elev = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
                /* Convert to radians: */
                double a = M_PI * (azim/180);
                double e = M_PI * (elev/180);
                /* Convert to unit direction vector: */
                src->dir.c[0] = cos(a)*cos(e);
                src->dir.c[1] = sin(a)*cos(e);
                src->dir.c[2] = sin(e);
                dir_given = TRUE;
              }
            else if (argparser_keyword_present_next(pp, "direction"))
              { if (dir_given) { argparser_error(pp, "lamp's direction" pst_bogus_spec_MESS); }
                /* Get Cartesian coordinates of direction vector: */
                r3_t vec;
                vec.c[0] = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
                vec.c[1] = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
                vec.c[2] = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
                if (r3_L_inf_norm(&vec) != 0.0)
                  { /* Convert to unit direction vector: */
                    (void)r3_dir(&vec, &(src->dir));
                  }
                dir_given = TRUE;
              }
            else if (argparser_keyword_present_next(pp, "radius"))
              { if (crad_given) { argparser_error(pp, "lamp's radius" pst_bogus_spec_MESS); }
                /* Get radius: */
                double radius = argparser_get_next_double(pp, 0.0, 180); 
                src->crad = cos(radius*M_PI/180);
                crad_given = TRUE;
              }
            else if (argparser_keyword_present_next(pp, "power"))
              { if (pwr_given) { argparser_error(pp, "lamp's power" pst_bogus_spec_MESS); }
                /* Get color: */
                src->pwr = frgb_parse_color(pp);
                pwr_given = TRUE;
              }
            else 
              { /* No prameter keyword next -- assume end of light spec: */
                break;
              }
          }
      }

    if (N_P != NULL) { (*N_P) = N; }
      
    return src;
  }

void pst_lamp_spec_write(FILE *wr, pst_lamp_t *src)
  { 
    fprintf(wr, "-lamp\n");
    r3_t *dir = &(src->dir);
    double crad = src->crad;
    frgb_t *pwr = &(src->pwr);
    
    if (r3_L_inf_norm(dir) != 0.0)
      { fprintf(wr, "  direction");
        fprintf(wr, " %+8.5f %+8.5f %+8.5f\n", dir->c[0], dir->c[1], dir->c[2]);
      }

    if (crad <= -1.0)
      { fprintf(wr, "  ambient\n"); }
    else if (crad == 0.0)
      { fprintf(wr, "  wall\n"); }
    else
      { fprintf(wr, "  radius %6.3f\n", acos(crad)); }
    
    if (! isnan(pwr->c[0]))
      { fprintf(wr, "  power");
        for (int32_t c = 0; c < 3; c++) { fprintf(wr, " %6.4f", pwr->c[c]); }
        fprintf(wr, "\n");
      }
    fflush(wr);
  }




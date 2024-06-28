/* See {boap_plate.h} */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <vec.h>
#include <jsfile.h>
#include <epswr.h>
#include <epswr_dim.h>

#include <boap_plate.h>

boap_plate_t *boap_plate_new(r3_t cmin, r3_t cmax, frgb_t *color)
  { boap_plate_t *plt = notnull(malloc(sizeof(boap_plate_t)), "no mem");
    plt->cmin = cmin;
    plt->cmax = cmax;
    plt->color = color;
    return plt;
  }

void boap_plate_ref_vec_draw 
  ( epswr_figure_t *epsf, 
    int32_t np,
    boap_plate_ref_vec_t *pltv, 
    bool_t fill, bool_t draw, 
    int8_t uax,
    double pos
  )
  {
    if (isnan(pos))
      { /* Render as 3D view with visibility. */
        /* Sort plates by increasing max {uax} coordinate: */
        /* !!! Should use a faster algo. !!! */
        for (int32_t i = 0; i < np-1; i++)
          { for (int32_t j = i+1; j < np; j++)
              { boap_plate_t *pti = pltv->e[i];
                boap_plate_t *ptj = pltv->e[j];
                if (pti->cmax.c[uax]> ptj->cmax.c[uax])
                  { /* Plate {j} has smaller {Xmax}, must be drawn before {i}: */
                    pltv->e[i] = ptj; 
                    pltv->e[j] = pti;
                  }
              }
          }
       /* Just draw plates: */
       for (int32_t i = 0; i < np; i++)  { boap_plate_draw(epsf, pltv->e[i], 0.25, fill, draw, uax, pos); }
      }
    else
      { /* Render as cross-sections: */
        if (draw)
          { for (int32_t i = 0; i < np; i++) { boap_plate_draw(epsf, pltv->e[i], 0.50, FALSE, TRUE, uax, pos); } }
        if (fill)
          { for (int32_t i = 0; i < np; i++) { boap_plate_draw(epsf, pltv->e[i], 0.00, TRUE, FALSE, uax, pos); } }
      }
  }
  
void boap_plate_draw 
  ( epswr_figure_t *epsf, 
    boap_plate_t *plt, 
    double penwd,
    bool_t fill, bool_t draw, 
    int8_t uax,
    double pos
  )
  {
    if ((! isnan(pos)) && ((pos < plt->cmin.c[uax]) || (pos > plt->cmax.c[uax]))) { return; }
    
    int8_t vax = (uax <= 0 ? 1 : 0);
    int8_t wax = (uax <= 1 ? 2 : 1);
    double xMin = plt->cmin.c[vax], xMax = plt->cmax.c[vax];
    double yMin = plt->cmin.c[wax], yMax = plt->cmax.c[wax];

    double R = (double)plt->color->c[0];
    double G = (double)plt->color->c[1];
    double B = (double)plt->color->c[2];
    epswr_set_fill_color(epsf, R,G,B);
    epswr_set_pen(epsf, 0,0,0, penwd, 0.0,0.0);
    epswr_rectangle(epsf, xMin,xMax, yMin,yMax, fill, draw);
  }

void boap_plate_print(FILE *wr, boap_plate_t *plt)
  {
    fprintf(wr, "  plate: ");
    fprintf(wr, " [ %7.1f _%7.1f]",   plt->cmin.c[0], plt->cmax.c[0]);
    fprintf(wr, " x [ %7.1f _%7.1f]", plt->cmin.c[1], plt->cmax.c[1]);
    fprintf(wr, " x [ %7.1f _%7.1f]", plt->cmin.c[2], plt->cmax.c[2]);
    fprintf(wr, "\n");
  }

vec_typeimpl(boap_plate_ref_vec_t,boap_plate_ref_vec,boap_plate_ref_t);

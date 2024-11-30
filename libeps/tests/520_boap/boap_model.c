/* See {boap_model.h} */

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
#include <jsfile.h>
#include <epswr.h>
#include <epswr_dim.h>

#include <boap_bar.h>
#include <boap_plate.h>

#include <boap_model.h>

boap_model_t *boap_model_new(void)
  {
    boap_model_t *model = notnull(malloc(sizeof(boap_model_t)), "no mem");
     
    model->nb = 0;
    model->barv = boap_bar_ref_vec_new(20);
   
    model->np = 0;
    model->pltv = boap_plate_ref_vec_new(20);
    
    model->nd = 0;
    model->dimv = boap_dim_ref_vec_new(20);

    return model;
  }

void boap_model_add_plate (boap_model_t *mod, boap_plate_t *plt)
  { 
    boap_plate_print(stderr, plt);
    boap_plate_ref_vec_expand(&(mod->pltv), mod->np);
    mod->pltv.e[mod->np] = plt;
    (mod->np)++;
    fprintf(stderr, "  np = %d\n", mod->np);
  }

void boap_model_add_bar (boap_model_t *mod, boap_bar_t *bar)
  { 
    boap_bar_print(stderr, bar);
    boap_bar_ref_vec_expand(&(mod->barv), mod->nb);
    mod->barv.e[mod->nb] = bar;
    (mod->nb)++;
    fprintf(stderr, "  nb = %d\n", mod->nb);
  }

void boap_model_add_dim (boap_model_t *mod, boap_dim_t *dim)
  { 
    boap_dim_print(stderr, dim);
    boap_dim_ref_vec_expand(&(mod->dimv), mod->nd);
    mod->dimv.e[mod->nd] = dim;
    (mod->nd)++;
  }

void boap_model_draw
  ( epswr_figure_t *epsf, 
    boap_model_t *mod, 
    int8_t uax,
    double pos
  )
  {
    /* Convert all bars to plates: */
    for (uint32_t ib = 0;  ib < mod->nb; ib++)
      { boap_bar_t *bar = mod->barv.e[ib];
        boap_bar_smash(bar, &(mod->pltv), &(mod->np));
      }
    fprintf(stderr, "  np = %d\n", mod->np);

    /* Draw the plates: */
    boap_plate_ref_vec_draw(epsf, mod->np, &(mod->pltv), TRUE, TRUE, uax, pos);
    
    /* Draw the dimension lines: */
    for (uint32_t id = 0;  id < mod->nd; id++)
      { boap_dim_t *dim = mod->dimv.e[id];
        boap_dim_draw(epsf, dim, uax);
      }
  }


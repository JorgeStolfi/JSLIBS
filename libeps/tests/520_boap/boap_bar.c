/* See {boap_bar.h} */

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

#include <boap_plate.h>

#include <boap_bar.h>

boap_bar_t *boap_bar_new
  ( char type,
    double Xctr,
    double Yctr,
    double Zctr,
    int8_t axA, double Asize, double Athk,
    int8_t axB, double Bsize, double Bthk,
    int8_t axC, double Csize,
    frgb_t *color
  )
  {
    demand ((axA >= 0) && (axA <= 2), "invalid A axis");
    demand ((axB >= 0) && (axB <= 2), "invalid B axis");
    demand ((axC >= 0) && (axC <= 2), "invalid C axis");
    demand ((axA != axB) && (axA != axC) && (axB != axC), "invalid axes");
    
    if (type != 'F')
      { demand (Asize*Athk > 0, "inconsistent signs {Asize,Athk}");
        demand (Bsize*Bthk > 0, "inconsistent signs {Bsize,Bthk}");
      }

    boap_bar_t *bar = notnull(malloc(sizeof(boap_bar_t)), "no mem");
    
    bar->type = type;
    
    bar->ctr = (r3_t){{ Xctr, Yctr, Zctr }};
    bar->wax = (i3_t){{ axA, axB, axC }};
    bar->size = (r3_t){{ Asize, Bsize, Csize }};
    bar->thk = (r2_t){{ Athk, Bthk }};
    bar->emin_trim = (r2_t){{ 0.0, 0.0 }};
    bar->emax_trim = (r2_t){{ 0.0, 0.0 }};
    bar->emin_acut = (r2_t){{ 0.0, 0.0 }};
    bar->emax_acut = (r2_t){{ 0.0, 0.0 }};
    
    /* !!! These defaults should depend on {Asize,Athk,Bsize,Bthk}: !!! */
    bar->irad = 1.0; /* Typical. */
    bar->orad = 0.5; /* Typical. */
    
    bar->color = color;
    
    return bar;
  }

void boap_bar_trim
  ( boap_bar_t *bar,
    int8_t bax,
    int8_t end,
    double amount
  )
  {
    demand ((bax >= 0) && (bax <= 1), "invalid part identifier {bax}");
    demand ((end == -1) || (end == +1), "invalid extremity identifier {end}");
    if (end < 0)
      { bar->emin_trim.c[bax] = amount;  }
    else
      { bar->emax_trim.c[bax] = amount;  }
  }

void boap_bar_smash
  ( boap_bar_t  *bar,
    boap_plate_ref_vec_t *plv,
    int32_t *npP
  )
  {
    int32_t np = (*npP);     /* Number of plates that make up the box. */
      
    double eps = 1.0e-4; /* Fudge displacement to avoid coincident surfaces. */

    auto void add_plate(int8_t bax, double Amin, double Amax, double Bmin, double Bmax);
      /* Adds to the plate set {plv} a part of {bar}, that extends 
        over {[Amin _ Amax]×[Bmin _ Bmax]} in the cross-section {A,B} coordinate system.
        Applies the trimming appropriate to the the {bax} parts of (0 for {A} parts, 
        or 1 for {B} parts). */
    
    void add_plate(int8_t bax, double Amin, double Amax, double Bmin, double Bmax)
      { 
        r3_t pmin, pmax; /* World coord ranges of plates: */
        
        /* World coord range along local {C} axis: */
        int8_t iwC = (int8_t)bar->wax.c[2]; /* World coord axis equated to local {C} coord axis. */
        double Csize = bar->size.c[2];
        pmin.c[iwC] = bar->ctr.c[iwC] - Csize/2 + bar->emin_trim.c[bax];
        pmax.c[iwC] = bar->ctr.c[iwC] + Csize/2 - bar->emax_trim.c[bax];
        
        /* World coord range along local {A} axis: */
        int8_t iwA = (int8_t)bar->wax.c[0]; /* World coord axis equated to local {A} coord axis. */
        pmin.c[iwA] = bar->ctr.c[iwA] + Amin;
        pmax.c[iwA] = bar->ctr.c[iwA] + Amax;
        
        /* World coord range along local {B} axis: */
        int8_t iwB = (int8_t)bar->wax.c[1]; /* World coord axis equated to local {B} coord axis. */
        pmin.c[iwB] = bar->ctr.c[iwB] + Bmin;
        pmax.c[iwB] = bar->ctr.c[iwB] + Bmax;
        
        /* Sort range ends: */
        for (int8_t k = 0; k < 3; k++)
          { if (pmax.c[k] < pmin.c[k]) 
              { double t = pmax.c[k]; pmax.c[k] = pmin.c[k]; pmin.c[k] = t; }
          }

        /* Shrink lightly: */
        for (int8_t k = 0; k < 3; k++) { pmin.c[k] += +eps; pmax.c[k] += -eps; }

        boap_plate_t *plt = boap_plate_new(pmin, pmax, bar->color);
        boap_plate_print(stderr, plt);
        boap_plate_ref_vec_expand(plv, np);
        plv->e[np] = plt;
        np++;
      }
        
    if (bar->type == 'F')
      { /* Flat bar -- origin is center of cross-section: */
        add_plate(0, -bar->size.c[0]/2, +bar->size.c[0]/2, -bar->thk.c[1]/2, +bar->thk.c[1]/2);
      }
    else if (bar->type == 'L')
      { /* 'L' bar -- origin is corner: */
        add_plate(0, 0, +bar->size.c[0], 0, bar->thk.c[1]);
        add_plate(1, 0, +bar->thk.c[0],  0, bar->size.c[1]);
      }
    else if (bar->type == 'T')
      { /* 'T' bar -- origin is top of 'T', but bar is upside down: */
        add_plate(0, -bar->size.c[0]/2, +bar->size.c[0]/2,  0,    +bar->thk.c[1] );
        add_plate(1, -bar->thk.c[0]/2,  +bar->thk.c[0]/2,   +eps, +bar->size.c[1]);
      }
    else if (bar->type == 'U')
      { /* 'U' bar -- origin is middle of bottom edge of base: */
        double Ard = bar->size.c[0]/2, Bsz = bar->size.c[1];
        double Ath = bar->thk.c[0], Bth = bar->thk.c[1];
        add_plate(0, -Ard+eps, +Ard-eps, 0,    Bth);
        add_plate(1, -Ard,     -Ard+Ath, +eps, +Bsz);
        add_plate(1, +Ard-Ath, +Ard,     +eps, +Bsz);
      }
    else if (bar->type == 'M')
      { /* 'M' bar -- origin is center of cross-section: */
        double Ard = bar->size.c[0]/2, Brd = bar->size.c[1]/2;
        double Ath = bar->thk.c[0], Bth = bar->thk.c[1];
        add_plate(0, -Ard+eps, +Ard-eps, -Brd,     -Brd+Bth);
        add_plate(0, -Ard+eps, +Ard-eps, +Brd-Bth, +Brd);
        add_plate(1, -Ard,     -Ard+Ath, -Brd+eps, +Brd-eps);
        add_plate(1, +Ard-Ath, +Ard,     -Brd+eps, +Brd-eps);
      }
    else
      { assert(FALSE); }
      
    (*npP) = np;
  }

void boap_bar_print(FILE *wr, boap_bar_t *bar)
  {
    fprintf(wr, "BAR %c", bar->type);
    fprintf(wr, " ctr = (%7.1f,%7.1f,%7.1f)", bar->ctr.c[0], bar->ctr.c[1], bar->ctr.c[2]);
    fprintf(wr, " axis %d: size = %7.2f thk = %7.2f", bar->wax.c[0], bar->size.c[0], bar->thk.c[0]);
    fprintf(wr, " axis %d: size = %7.2f thk = %7.2f", bar->wax.c[1], bar->size.c[1], bar->thk.c[1]);
    fprintf(wr, " axis %d: size = %7.2f", bar->wax.c[2], bar->size.c[2]);
    fprintf(wr, "\n");
  }

vec_typeimpl(boap_bar_ref_vec_t, boap_bar_ref_vec, boap_bar_ref_t);

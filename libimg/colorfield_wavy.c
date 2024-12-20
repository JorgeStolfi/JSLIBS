/* see colorfield_wavy.h */
/* Last edited on 2024-12-04 22:48:09 by stolfi */

#define colorfield_wavy_COPYRIGHT "(C) 2003 by Jorge Stolfi, the University of Campinas, Brazil."

/* See the rights and conditions notice at the end of this file. */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <jsmath.h>
#include <bool.h>
#include <argparser.h>

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>
#include <colorfield_wavy.h>

/* IMPLEMENTATIONS */

#define TwoPi (2.0 * M_PI)

cfld_wavy_args_t *cfld_wavy_parse_uniform(argparser_t *pp) 
  {
    cfld_wavy_args_t *wfa = (cfld_wavy_args_t *)malloc(sizeof(cfld_wavy_args_t));
    frgb_t color = frgb_parse_color(pp); 
    wfa->org = (cfld_int_pair_t){{0,0}};
    wfa->orgColor = color;
    wfa->waves = NULL;
    return wfa;
  }

cfld_wavy_args_t *cfld_wavy_parse_simple(argparser_t *pp, int32_t uX, int32_t uY) 
  {
    cfld_wavy_args_t *wfa = (cfld_wavy_args_t *)malloc(sizeof(cfld_wavy_args_t));
    int32_t phase0 = (int32_t)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
    frgb_t color0 = frgb_parse_color(pp); 
    int32_t phase1 = (int32_t)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
    frgb_t color1 = frgb_parse_color(pp);
    cfld_wave_args_t *wa = (cfld_wave_args_t *)malloc(sizeof(cfld_wave_args_t));

    /* Max and min phases must be different: */
    if(phase0 == phase1)
      { argparser_error(pp, "max and min positions must be different"); }

    /* Fill in the wave's descriptor: */
    wa->bot = (cfld_int_pair_t){{phase1*uX, phase1*uY}};
    wa->botColor = color1;
    wa->next = NULL;

    /* Fill in the wavy field descriptor: */
    wfa->org = (cfld_int_pair_t){{phase0*uX, phase0*uY}};
    wfa->orgColor = color0;
    wfa->waves = wa;
    
    return wfa;
  }
  
cfld_wavy_args_t *cfld_wavy_parse_general(argparser_t *pp, int32_t nWaves)
  {
    cfld_wavy_args_t *wfa = (cfld_wavy_args_t *)malloc(sizeof(cfld_wavy_args_t));
    int32_t X0 = (int32_t)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
    int32_t Y0 = (int32_t)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
    frgb_t color0 = frgb_parse_color(pp);
    
    wfa->orgColor = color0;
    wfa->org = (cfld_int_pair_t){{X0, Y0}};
    wfa->waves = NULL;
    for (int32_t k = 0; k < nWaves; k++)
      { cfld_wave_args_t *wa = (cfld_wave_args_t *)malloc(sizeof(cfld_wave_args_t));
        int32_t Xk = (int32_t)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
        int32_t Yk = (int32_t)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
        frgb_t colork = frgb_parse_color(pp);
        
        /* Period must be nonzero: */
        if((Xk == X0) && (Yk == Y0))
          { argparser_error(pp, "max and min positions must be different"); }
        
        wa->bot = (cfld_int_pair_t){{Xk, Yk}};
        wa->botColor = colork; 
        wa->next = wfa->waves;
        wfa->waves = wa;
      }
      
    return wfa;
  }

cfld_wavy_params_t *cfld_wavy_compute_params
  ( cfld_wavy_args_t *wfa, 
    frgb_adjuster_t *adjust,
    bool_t logarithmic
  )
  { cfld_wavy_params_t *wfp = (cfld_wavy_params_t *)malloc(sizeof(cfld_wavy_params_t));
    /* Note: user-input colors are all RGB, even when {o->gray} is true. */
    
    wfp->org = wfa->org;
    wfp->orgColor = adjust(&(wfa->orgColor), wfa->org.c[0], wfa->org.c[1]);
    frgb_print(stderr, "orgColor = ( ", &(wfa->orgColor), 3, "%7.5f", " )\n");
    
    /* Compute wave-adjustment tables: */
    { cfld_wave_args_t *wa = wfa->waves;
      wfp->waves = NULL;
      while (wa != NULL)
        { frgb_t botColor = adjust(&(wa->botColor), wa->bot.c[0], wa->bot.c[1]);
          cfld_wave_params_t *wp = cfld_wavy_compute_wave_params
            ( wa, &(wfp->org), &(wfp->orgColor), &botColor, logarithmic );
          wp->next = wfp->waves;
          wfp->waves = wp;
          wa = wa->next;
        }
    }
    return wfp;
  }
 
cfld_wave_params_t *cfld_wavy_compute_wave_params
  ( cfld_wave_args_t *wa,
    cfld_int_pair_t *org,
    frgb_t *orgColor, 
    frgb_t *botColor,
    bool_t logarithmic
  )
  {
    cfld_wave_params_t *wavep = (cfld_wave_params_t *)malloc(sizeof(cfld_wave_params_t));
    /* Compute period vector {dX,dY}: */
    int32_t dX = 2*(wa->bot.c[0] - org->c[0]);
    int32_t dY = 2*(wa->bot.c[1] - org->c[1]);
    int32_t g = (int32_t)gcd((uint64_t)abs(dX), (uint64_t)abs(dY));
    cfld_int_pair_t frN = (cfld_int_pair_t){{dX/g, dY/g}};
    int32_t frD = dX*frN.c[0] + dY*frN.c[1];
    wavep->freqNum = frN;
    wavep->freqDen = frD;
    
    { frgb_t *tb = talloc(frD, frgb_t);
      int32_t iphase, i;
      frgb_t ampl;
      
      /* Debugging: */
      fprintf(stderr, "\n");
      frgb_print(stderr, "orgColor = ( ", orgColor, 3, "%7.5f", " )\n");
      frgb_print(stderr, "botColor = ( ", botColor, 3, "%7.5f", " )\n");
      fprintf(stderr, "logarithmic = %c\n", "FT"[logarithmic]);
      
      /* Compute the peak-to-peak amplitude of the wave: */
      for (i = 0; i < 3; i++) 
        { if (logarithmic) 
            { double boti = botColor->c[i];
              double orgi = orgColor->c[i];
              if (boti < VAL_EPS) { boti = VAL_EPS; }
              if (orgi < VAL_EPS) { orgi = VAL_EPS; }
              ampl.c[i] = (float)log(boti / orgi);
            }
         else
           { ampl.c[i] = (float)(botColor->c[i] - orgColor->c[i]); }
        }
      /* Compute wave values at each {iphase}: */
      for (iphase = 0; iphase < frD; iphase++)
        {
          frgb_t tbi; 
          double phase = ((double)iphase)/frD;
          double wfun = 0.5 * (1.0 - cos(TwoPi*phase));
          for (i = 0; i < 3; i++) 
            { if (logarithmic) 
                { tbi.c[i] = (float)exp(ampl.c[i] * wfun); }
              else
                { tbi.c[i] = (float)(ampl.c[i] * wfun); }
            }
          /* Debugging: */
          if ((iphase < 5 ) || (iphase % (frD/10) == 0) || (iphase >= frD-5))
            { fprintf(stderr, "tb[%05d] = ", iphase);
              frgb_print(stderr, "( ", &tbi, 3, "%7.5f", " )\n");
            }
          
          tb[iphase] = tbi;
        }
      wavep->tb = tb;
    }
    wavep->next = NULL;
    return wavep;
  }
  
void cfld_wavy_eval
  ( cfld_wavy_params_t *wfp,
    bool_t logarithmic, 
    int32_t col, 
    int32_t row,
    frgb_t *fv,
    int32_t chns
  )
  {
    *fv = wfp->orgColor;
    if (wfp->waves != NULL)
      { int32_t dCol = col - wfp->org.c[0];
        int32_t dRow = row - wfp->org.c[1];
        cfld_wavy_apply_waves(fv, chns, dCol, dRow, wfp->waves, logarithmic);
      }
  }

void cfld_wavy_apply_waves
  ( frgb_t *locColor,
    int32_t chns,
    int32_t dCol, 
    int32_t dRow,
    cfld_wave_params_list_t wp,
    bool_t logarithmic
  )
  {
    int32_t i;
    while (wp != NULL)
      { frgb_t *tb = wp->tb;
        cfld_int_pair_t frN = wp->freqNum;
        int32_t iphase = (int32_t)imod((dCol*frN.c[0] + dRow*frN.c[1]), wp->freqDen);
        float *tbxy = tb[iphase].c;
        for (i = 0; i < chns; i++)
          { if (logarithmic) 
              { locColor->c[i] *= tbxy[i]; }
            else
              { locColor->c[i] += tbxy[i]; }
          }
        wp = wp->next;
      }
  }

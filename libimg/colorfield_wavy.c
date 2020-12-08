/* see colorfield_wavy.h 
** Last edited on 2017-01-02 16:27:45 by jstolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <jsmath.h>
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

cfld_wavy_args_t *cfld_wavy_parse_simple(argparser_t *pp, int uX, int uY) 
  {
    cfld_wavy_args_t *wfa = (cfld_wavy_args_t *)malloc(sizeof(cfld_wavy_args_t));
    int phase0 = (int)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
    frgb_t color0 = frgb_parse_color(pp); 
    int phase1 = (int)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
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
  
cfld_wavy_args_t *cfld_wavy_parse_general(argparser_t *pp, int nWaves)
  {
    cfld_wavy_args_t *wfa = (cfld_wavy_args_t *)malloc(sizeof(cfld_wavy_args_t));
    int k;
    int X0 = (int)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
    int Y0 = (int)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
    frgb_t color0 = frgb_parse_color(pp);
    
    wfa->orgColor = color0;
    wfa->org = (cfld_int_pair_t){{X0, Y0}};
    wfa->waves = NULL;
    for (k = 0; k < nWaves; k++)
      { cfld_wave_args_t *wa = (cfld_wave_args_t *)malloc(sizeof(cfld_wave_args_t));
        int Xk = (int)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
        int Yk = (int)argparser_get_next_int(pp, -INT_MAX, INT_MAX);
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
    int logarithmic
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
    int logarithmic
  )
  {
    cfld_wave_params_t *wavep = (cfld_wave_params_t *)malloc(sizeof(cfld_wave_params_t));
    /* Compute period vector {dX,dY}: */
    int dX = 2*(wa->bot.c[0] - org->c[0]);
    int dY = 2*(wa->bot.c[1] - org->c[1]);
    int g = (int)gcd(abs(dX), abs(dY));
    cfld_int_pair_t frN = (cfld_int_pair_t){{dX/g, dY/g}};
    int frD = dX*frN.c[0] + dY*frN.c[1];
    wavep->freqNum = frN;
    wavep->freqDen = frD;
    
    { frgb_t *tb = (frgb_t *)malloc(frD * sizeof(frgb_t));
      int iphase, i;
      frgb_t ampl;
      
      /* Debugging: */
      fprintf(stderr, "\n");
      frgb_print(stderr, "orgColor = ( ", orgColor, 3, "%7.5f", " )\n");
      frgb_print(stderr, "botColor = ( ", botColor, 3, "%7.5f", " )\n");
      fprintf(stderr, "logarithmic = %d\n", logarithmic);
      
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
    int logarithmic, 
    int col, 
    int row,
    frgb_t *fv,
    int chns
  )
  {
    *fv = wfp->orgColor;
    if (wfp->waves != NULL)
      { int dCol = col - wfp->org.c[0];
        int dRow = row - wfp->org.c[1];
        cfld_wavy_apply_waves(fv, chns, dCol, dRow, wfp->waves, logarithmic);
      }
  }

void cfld_wavy_apply_waves
  ( frgb_t *locColor,
    int chns,
    int dCol, 
    int dRow,
    cfld_wave_params_list_t wp,
    int logarithmic
  )
  {
    int i;
    while (wp != NULL)
      { frgb_t *tb = wp->tb;
        cfld_int_pair_t frN = wp->freqNum;
        int iphase = (int)imod((dCol*frN.c[0] + dRow*frN.c[1]), wp->freqDen);
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

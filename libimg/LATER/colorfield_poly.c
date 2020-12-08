/* See colorfield_poly.h 
** Last edited on 2008-05-25 03:22:26 by stolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

/* TO BE COMPLETED */

#include <colorfield_poly.h>

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include <argparser.h>
#include <jsmath.h>

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>

cfld_poly_params_t *cfld_poly_compute_bezier_coeffs
  ( cfld_poly_args_t *w, 
    frgb_t *orgColor, 
    frgb_t *botColor
  )
  {
    cfld_poly_params_t *wtb_el = (cfld_poly_params_t *)malloc(sizeof(cfld_poly_params_t));
    int g = gcd(abs(w->period.c[0]), abs(w->period.c[1]));
    cfld_int_pair_t frN = (cfld_int_pair_t){{w->period.c[0]/g, w->period.c[1]/g}};
    int frD = w->period.c[0]*frN.c[0] + w->period.c[1]*frN.c[1];
    wtb_el->freqNum = frN;
    wtb_el->freqDen = frD;
    wtb_el->additive = w->additive;
    
    { frgb_t *tb = (frgb_t *)malloc(frD * sizeof(frgb_t));
      int iphase, i;
      frgb_t ampl;
      fprintf(stderr, "\n");
      print_triplet(stderr, "orgColor = ( ", 3, orgColor, " )\n");
      print_triplet(stderr, "botColor = ( ", 3, botColor, " )\n");
      fprintf(stderr, "additive = %d\n", w->additive);
      
      /* Compute the peak-to-peak amplitude of the wave: */
      for (i = 0; i < 3; i++) 
        { if (w->additive) 
            { ampl.c[i] = botColor->c[i] - orgColor->c[i]; }
          else
            { ampl.c[i] = log(botColor->c[i] / orgColor->c[i]); }
        }
      /* Compute wave values at each {iphase}: */
      for (iphase = 0; iphase < frD; iphase++)
        {
          frgb_t tbi; 
          double phase = ((double)iphase)/frD;
          double wfun = 0.5 * (1.0 - cos(TwoPi*phase));
          for (i = 0; i < 3; i++) 
            { if (w->additive) 
                { tbi.c[i] = ampl.c[i] * wfun; }
              else
                { tbi.c[i] = exp(ampl.c[i] * wfun); }
            }
          if ((iphase < 5 ) || (iphase % (frD/10) == 0) || (iphase >= frD-5))
            { fprintf(stderr, "tb[%05d] = ", iphase);
              print_triplet(stderr, "( ", 3, &tbi, " )\n");
            }
          tb[iphase] = tbi;
        }
      wtb_el->tb = tb;
    }
    wtb_el->next = NULL;
    return wtb_el;
  }
  
void cfld_poly_eval_bezier
  ( frgb_t *locColor,
    cfld_poly_params_t wtb,
    int dCol, 
    int dRow
  )
  {
    int i;
    while (wtb != NULL)
      { frgb_t *tb = wtb->tb;
        cfld_int_pair_t frN = wtb->freqNum;
        int iphase = imod((dCol*frN.c[0] + dRow*frN.c[1]), wtb->freqDen);
        frgb_t tbxy = tb[iphase];
        for (i = 0; i < 3; i++)
          { if (wtb->additive) 
              { locColor->c[i] += tbxy.c[i]; }
            else
              { locColor->c[i] *= tbxy.c[i]; }
          }
        wtb = wtb->next;
      }
  }

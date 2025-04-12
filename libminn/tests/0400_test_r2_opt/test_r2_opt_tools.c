/* See {test_r2_opt_tools.h}. */
/* Last edited on 2025-03-19 14:19:04 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <affirm.h>

#include <test_r2_opt_basic.h>

#include <test_r2_opt_tools.h>

double tr2o_tools_eval_mother_image
  ( r2_t *p, 
    r2_t *scale,
    uint32_t mom_NF, 
    r2_t mom_frq[], 
    r2_t mom_phi[], 
    double mom_amp[]
  )
  {
    int32_t t;
    double x = p->c[0];      /* Point (unscaled) in X. */
    double y = p->c[1];      /* Point (unscaled) in Y. */
    double sx = scale->c[0]; /* Image scaling factor in X. */
    double sy = scale->c[1]; /* Image scaling factor in Y. */
    double rf0 = 0.25;       /* Cutoff frequency (cycles per shrunk pixel). */
    double sum_val = 0.0;    /* Sum ov values of all component waves. */
    double sum_amp = 0.0;    /* Sum of amplitude of all component waves. */
    for (t = 0; t < 10; t++)
      { /* Get the frequency vector of term {t} in shrunken image: */
        double fx = mom_frq[t].c[0];   /* X frequency (cycles per original pixel). */
        double fy = mom_frq[t].c[1];   /* Y frequency (cycles per original pixel). */
        /* Wave phases of each factor of term {t} at {x,y} (as fraction of cycle). */
        double tph0 = (+ fx*x + fy*y) + mom_phi[t].c[0]; 
        double tph1 = (- fy*x + fx*y) + mom_phi[t].c[1]; 
        /* Raw vale of term {t}: */
        double raw = sin(2*M_PI*tph0)*sin(2*M_PI*tph1); 
        /* Frequency in cycles per srunk pixel: */
        double rfx = fx*sx;
        double rfy = fy*sy;
        /* Attenuation factor of term {t}: */
        double fr2 = (rfx*rfx + rfy*rfy)/(rf0*rf0); /* Freq rel to cutoff, squared. */
        double att = (fr2 >= 1.0 ? 0.0 : 1 - sin(0.5*M_PI*fr2)); /* Attenuation factor. */
        /* fprintf(stderr, "!! %3d fx = %12.8f  fy = %12.8f  fr2 = %12.8f  att = %12.8f\n", t, fx, fy, fr2, att); */
        assert((att >= 0) && (att <= 1.0));
        /* Amplitude of term {t} after attenuation: */
        double amp = att*mom_amp[t]; 
        assert(amp >= 0);
        /* Value of term {t}: */
        double val = amp*raw;
        /* Accumulate: */
        sum_val += val;
        sum_amp += amp;
      }
    /* Map to range {[0 _ 1]}: */
    /* fprintf(stderr, "!! %12.8f %12.8f\n", sum_val, sum_amp); */
    /* assert(fabs(sum_val) <= sum_amp); */
    return (sum_val/sum_amp + 1)/2;
  }
      
void tr2o_tools_initialize_mother_image
  ( uint32_t mom_NF, 
    r2_t mom_frq[], 
    r2_t mom_phi[], 
    double mom_amp[],
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "  mother function terms:\n"); }
    demand(mom_NF >= 2, "need at least 2 waves");
    double gold = 0.6180398876; /* Golden ratio. */
    double Pmin =  2.5;       /* Min wavelength (pixels). */
    double Pmax = 1024*Pmin;  /* Max wavelength (pixels). */
    double Pstep = exp(log(Pmax/Pmin)/(mom_NF-1));
    for (uint32_t t = 0;  t < mom_NF; t++)
      { /* Choose the frequency vector of term {t}: */
        double P = Pmin*pow(Pstep,t);     /* Wavelength at scale 0 (pixels). */
        double azm = gold*2*M_PI*t;       /* Azimuth of wave direction (radians). */ 
        double fx = cos(azm)/P;           /* X frequency at scale 0 (cycles per pixel). */
        double fy = sin(azm)/P;           /* Y frequency at scale 0 (cycles per pixel). */
        mom_frq[t] = (r2_t){{ fx, fy }};  /* Frequency vector at scale 0 (cycles per pixel). */
        for (uint32_t r = 0;  r < 2; r++)
          { double phi = gold*(2*t+7+r);          /* Wave initial phase (as fraction of cycle). */
            mom_phi[t].c[r] = phi - floor(phi);   /* Reduce to {[0 _ 1]}. */
          }
        mom_amp[t] = 1.0;                 /* Raw amplitude. 1.5 + cos(3.5*t); */
        if (verbose)
          { fprintf(stderr, "    %3d  frq = ( %10.5f %10.5f )", t,  mom_frq[t].c[0], mom_frq[t].c[1]);
            fprintf(stderr, "  per = %10.5f", P);
            fprintf(stderr, "  phi = ( %10.5f %10.5f ) amp = %10.5f\n", mom_phi[t].c[0], mom_phi[t].c[1], mom_amp[t]);
          }
      }
  }

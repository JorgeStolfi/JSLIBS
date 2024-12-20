/* See epswr_dim.h */
/* Last edited on 2024-12-05 10:13:45 by stolfi */

#include <math.h>

#include <affirm.h>

#include <epswr.h>
#include <epswr_dev.h>
#include <epswr_def.h>
#include <epswr_dev_dim.h>

void epswr_dev_dim_linear
  ( epswr_figure_t *eps, 
    double psxa, double psya,
    double psxb, double psyb,
    double psagap,
    double psbgap, 
    double pselen,
    bool_t inner,
    double psdpos, double psdlen, 
    double pshoff, double psvoff,
    double *psxrP, double *psyrP,
    double *rotP
  )
  { double psxd = psxb - psxa;
    double psyd = psyb - psya;
    double psdab = hypot(psxd, psyd);
    demand(psdab > 0.0, "coincident points");
    
    /* Compute unit vector {u=(psxu,psyu)} parallel to {a--b}: */
    double psxu = psxd/psdab;
    double psyu = psyd/psdab;
    
    /* Compute unit vector {v=(psxv,psyv)} at +90 degrees from {u}: */
    double psxv = -psyu;
    double psyv = psxu;
    
    if (eps->verbose)
      { fprintf(stderr, "psa = (%.2f,%.2f) psb =  (%.2f,%.2f)\n", psxa, psya, psxb, psyb); }
    
    /* Extremal signed distance of extension start in the direction of {pselen}: */
    double psmgap = (pselen > 0.0 ?  fmax(psagap, psbgap) : fmin(psagap, psbgap));

    /* Signed distance from extension ends to line {a--b}: */
    double psetip = psmgap + pselen; 
    
    /* Signed distance from dimension segment to line {a--b}: */
    double psdgap = (pselen > 0.0 ? psetip - psdpos : psetip + psdpos);

    if (pselen != 0)
      { 
        /* Extension segment on {A}: */
        double psxa0 = psxa + psagap*psxv;
        double psya0 = psya + psagap*psyv;
        double psxa1 = psxa + psetip*psxv;
        double psya1 = psya + psetip*psyv;
        epswr_dev_segment(eps, psxa0, psya0, psxa1, psya1);

        /* Extension segment on {B}: */
        double psxb0 = psxb + psbgap*psxv;
        double psyb0 = psyb + psbgap*psyv;
        double psxb1 = psxb + psetip*psxv;
        double psyb1 = psyb + psetip*psyv;
        epswr_dev_segment(eps, psxb0, psyb0, psxb1, psyb1);
      }
    
    if (psdlen > 0)
      { /* Draw the dimension segment arrow(s). */
        
        /* Tips of the arrows: */
        double psxa1 = psxa + psdgap*psxv;
        double psya1 = psya + psdgap*psyv;
        double psxb1 = psxb + psdgap*psxv;
        double psyb1 = psyb + psdgap*psyv;

        /* Bases of the arrows: */
        double psxa0, psya0, psxb0, psyb0;

        if (inner)
          { /* Draw arrow stems between line {A.B}: */ 
            if (psdlen >= psdab/2)
              { /* Single double-headed arrow: */
                psxa0 = psxb1;
                psya0 = psyb1;
                psxb0 = psxa1;
                psyb0 = psya1;
                epswr_dev_segment(eps, psxa1, psya1, psxb1, psyb1);
              }
            else
              { /* Two outward-pointing arrows: */
                psxa0 = psxa1 + psdlen*psxu;
                psya0 = psya1 + psdlen*psyu;
                psxb0 = psxb1 - psdlen*psxu;
                psyb0 = psyb1 - psdlen*psyu;
                epswr_dev_segment(eps, psxa0, psya0, psxa1, psya1);
                epswr_dev_segment(eps, psxb0, psyb0, psxb1, psyb1);
              }
          }
        else
          { /* Two inward-pointing arrows: */
            psxa0 = psxa1 - psdlen*psxu;
            psya0 = psya1 - psdlen*psyu;
            psxb0 = psxb1 + psdlen*psxu;
            psyb0 = psyb1 + psdlen*psyu;
            epswr_dev_segment(eps, psxa0, psya0, psxa1, psya1);
            epswr_dev_segment(eps, psxb0, psyb0, psxb1, psyb1);
          }
        /* Draw arrowheads: */
        double pslength = 7.5; /* Arrow length (pt). */
        double pswidth = 2.0;  /* Arrow width (pt). */
        epswr_dev_arrowhead
          ( eps, psxa0, psya0, psxa1, psya1,
            pswidth/2, pswidth/2, pslength, 1.0, 
            TRUE, TRUE
          );
        epswr_dev_arrowhead
          ( eps, psxb0, psyb0, psxb1, psyb1,
            pswidth/2, pswidth/2, pslength, 1.0, 
            TRUE, TRUE
          );
      }
    /* Label reference point: */
    if (psxrP != NULL) { (*psxrP) = (psxa + psxb)/2 + (psdgap + psvoff)*psxv + pshoff*psxu; }
    if (psyrP != NULL) { (*psyrP) = (psya + psyb)/2 + (psdgap + psvoff)*psyv + pshoff*psyu; }
    if (rotP != NULL) { (*rotP) = 180/M_PI*atan2(psyd, psxd); }
  }    
    

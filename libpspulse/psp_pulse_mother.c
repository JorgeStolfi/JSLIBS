/* See psp_pulse_mother.h */
/* Last edited on 2011-09-19 01:04:18 by stolfilocal */

#include <math.h>

#include <affirm.h>
#include <bz_basic.h>

#include <psp_basic.h>
#include <psp_pulse_mother.h>
#include <psp_pulse_mother_B.h>
#include <psp_pulse_mother_H.h>
#include <psp_pulse_mother_N.h>

psp_pulse_mother_count_t psp_pulse_mother_count
  ( psp_pulse_kind_t pkind, 
    psp_cont_t c, 
    psp_degree_t g
  )
  { unsigned int nmp;
    switch(pkind)
      { case psp_PK_B: 
          nmp = psp_pulse_mother_B_count(c, g); break;
        case psp_PK_H: 
          nmp = psp_pulse_mother_H_count(c, g); break;
        case psp_PK_N: 
          nmp = psp_pulse_mother_N_count(c, g); break;
        default:
          affirm(FALSE, "unknown pulse kind");
          nmp = 0;
      }
    return nmp; 
  }

psp_grid_size_t psp_pulse_mother_supp_count
  ( psp_pulse_family_t *fam, 
    psp_pulse_mother_index_t pix
  )
  { switch(fam->pkind)
      { case psp_PK_B: 
          return psp_pulse_mother_B_supp_count(fam->c, fam->g, pix);
        case psp_PK_H: 
          return psp_pulse_mother_H_supp_count(fam->c, fam->g, pix);
        case psp_PK_N: 
          return psp_pulse_mother_N_supp_count(fam->c, fam->g, pix);
        default:
          affirm(FALSE, "unknown pulse kind");
          return 0.0; /* To keep the compiler happy. */
      }
  }

void psp_pulse_mother_eval
  ( psp_pulse_family_t *fam, 
    psp_pulse_mother_index_t pix, 
    double z,
    psp_cont_t ord,
    double *w
  )
  { affirm(fam->g >= 0, "negative degree");
    affirm(fam->g <= psp_pulse_MAX_DEGREE, "excessive degree");
    affirm(ord <= psp_pulse_MAX_DIFF_ORD, "diff order too high");
    /* Get B�zier coeffs for interval containing {z}: */
    double bz[psp_pulse_MAX_DEGREE + 1];
    int k = (int) floor(z); 
    psp_pulse_mother_to_bezier(fam, pix, k, bz);
    /* Apply DeCasteljau's algorithm: */
    bz_eval(fam->g, bz, 0.0, 1.0, z - k, ord, w);
  }

void psp_pulse_mother_to_bezier
  ( psp_pulse_family_t *fam, 
    psp_pulse_mother_index_t pix,
    psp_grid_pos_t ix,
    double bz[]
  )
  { affirm(fam->g >= 0, "negative degree");
    affirm(fam->g <= psp_pulse_MAX_DEGREE, "excessive degree");
    switch(fam->pkind)
      { case psp_PK_B: 
          psp_pulse_mother_B_to_bezier(fam->c, fam->g, pix, ix, bz);
          return;
        case psp_PK_H: 
          psp_pulse_mother_H_to_bezier(fam->c, fam->g, pix, ix, bz);
          return;
        case psp_PK_N: 
          psp_pulse_mother_N_to_bezier(fam->c, fam->g, pix, ix, bz);
          return;
        default:
          affirm(FALSE, "unknown pulse kind");
      }
  }

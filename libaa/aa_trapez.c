/* See {aa_trapez.h}. */
/* Last edited on 2008-01-18 21:41:32 by stolfi */

#include <aa_trapez.h>

#include <ia.h>
#include <ia_trapez.h>
#include <aa.h>
#include <affirm.h>

ia_trapez_t aa_trapez_from_pair(Interval *xr, AAP xf, AAP yf)
  { ia_trapez_t tr;
    tr.x = *xr;
    if (xf->nterms == 0)
      { /* {xf} is a constant. */
        tr.yxlo = aa_range(yf);
        tr.yxhi = tr.yxlo;
      }
    else
      { demand((xf->nterms == 1), "xf must depend on only one noise variable");
        /* Get the midpoint of {xf}'s implicit range: */
        Float xmid = xf->center;

        /* Get index {id} and coef {xrad} of {xf}'s noise variable: */
        VarId id = ((AATermP)(xf + 1))->id;
        Float xrad = ((AATermP)(xf + 1))->coef;

        /* Get coef {yrad} of noise var {id} in {yf}: */
        AATerm eps[1];
        eps[0].id = id;
        aa_get_eps_coef(yf, 1, eps);
        Float yrad = eps[0].coef;
        
        /* Compute Y range {yxmid} at {xmid}: */
        MemP frame = aa_top();
        eps[0].coef = Zero;
        Interval yxmid = aa_implicit_range(aa_fix_eps(yf, 1, eps));
        aa_flush(frame);
        
        /* Compute slope interval {dyr} for midline: */
        Interval dyr = ia_div(ia_const(yrad,0), ia_const(xrad,0));

        /* Extrapolate the {Y} range from {xmid} to the ends of {xr}: */
        Interval dxlo = ia_sub(ia_const(xr->lo,0), ia_const(xmid,0));
        tr.yxlo = ia_add(yxmid, ia_mul(dxlo, dyr));
        
        Interval dxhi = ia_sub(ia_const(xr->hi,0), ia_const(xmid,0));
        tr.yxhi = ia_add(yxmid, ia_mul(dxhi, dyr));
      }
      
    return tr;
  }

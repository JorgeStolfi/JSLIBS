/* See cpk_weight.h */
/* Last edited on 2024-12-31 16:37:32 by stolfi */

#include <stdint.h>
#include <math.h>

#include <cpk_weight.h>

/* Range of each weight term after quantization. */
#define MZ ((double)10)

/* Weight multipliers for tiebreaking terms, by priority */
#define M1 ((double)10000)
#define M2 ((double)10000)
#define M3 ((double)1000)

/* Normalization factor so that priority 0 terms are in natural units: */
#define MM (MZ*M1*M2*M3)

double cpk_compute_weight(cpk_wcoeffs_t *wc, double URB, double DEM, double CLS)
  {
    double Z[4];
    Z[0] = 0; Z[1] = 0; Z[2] = 0; Z[3] = 0;
    if (wc->pNUM < 4) { Z[wc->pNUM] += wc->cNUM; }
    if (wc->pURB < 4) { Z[wc->pURB] += wc->cURB*URB; }
    if (wc->pDEM < 4) { Z[wc->pDEM] += wc->cDEM*DEM; }
    if (wc->pCLS < 4) { Z[wc->pCLS] += wc->cCLS*CLS; }

    /*  The weight {W} is obtained by computing {K[p] = round(10*Z[p])/10}
    and taking {W = K[0] + K[1]/10^4 + K[2]/10^8 + K[3]/10^11}. This
    trick will give the expected ranking of solutions as long as {nJ} is
    less than 1000. For {nJ} between 1000 and 10000, a large difference
    in {K[3]} may overcome small contrary differences in {K[2]}, and the
    lower bits of {K[3]} may be lost due to roundoff. For {nJ} above
    10000, the same reversal may occur at the other levels. */

    { int32_t p; for (p = 0; p < 4; p++) { Z[p] = floor(MZ*Z[p] + 0.5); } }
    
    return (((Z[0]*M1 + Z[1])*M2 + Z[2])*M3 + Z[3])/MM;
  }

void cpk_normalize_weight_coeffs(cpk_wcoeffs_t *wc)
  { /* Compute coefficient sum {S[p]} by priority {p}: */
    double S[4];
    S[0] = 0; S[1] = 0; S[2] = 0; S[3] = 0;
    if (wc->pNUM < 4) { S[wc->pNUM] += wc->cNUM; }
    if (wc->pURB < 4) { S[wc->pURB] += wc->cURB; }
    if (wc->pDEM < 4) { S[wc->pDEM] += wc->cDEM; }
    if (wc->pCLS < 4) { S[wc->pCLS] += wc->cCLS; }
    /* Normalize coeffs, set irrelevent ones to 0: */
    if (wc->pNUM < 4) { wc->cNUM /= S[wc->pNUM]; } else { wc->cNUM = 0; }
    if (wc->pURB < 4) { wc->cURB /= S[wc->pURB]; } else { wc->cURB = 0; }
    if (wc->pDEM < 4) { wc->cDEM /= S[wc->pDEM]; } else { wc->cDEM = 0; }
    if (wc->pCLS < 4) { wc->cCLS /= S[wc->pCLS]; } else { wc->cCLS = 0; }
  }


void cpk_print_weights(FILE *wr, cpk_wcoeffs_t *wc)
  { 
    auto void prt(char *name, uint32_t p, double c);
  
    for (uint32_t p = 0; p < 4; p++) 
      { if (wc->pNUM == p) { prt("NUM", p, wc->cNUM); }
        if (wc->pURB == p) { prt("URB", p, wc->cURB); }
        if (wc->pDEM == p) { prt("DEM", p, wc->cDEM); }
        if (wc->pCLS == p) { prt("CLS", p, wc->cCLS); }
      }
      
    void prt(char *name, uint32_t p, double c)
      { fprintf(wr, "  %s priority %d coeff = %8.6f\n", name, p, c); }
  }

int32_t cpk_weight_cmp(double Wa, double Wb)
  { double tol = 1.0e-14*(fabs(Wa)/2 + fabs(Wb)/2);
    if (fabs(Wa - Wb) <= tol)
      { return 0; }
    else if (Wa < Wb)
      { return -1; }
    else
      { return +1; }
  }

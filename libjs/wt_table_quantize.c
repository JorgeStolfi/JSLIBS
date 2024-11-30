/* See wt_table_quantize.h */
/* Last edited on 2024-11-19 04:39:00 by stolfi */

#define wt_table_quantize_C_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <jsmath.h>

#include <wt_table_quantize.h>

uint64_t wt_table_quantize
  ( uint32_t n,
    double wf[],
    uint64_t wia_sum,
    bool_t keep_nz,
    int32_t wi[]
  )
  { 
    bool_t debug = (n < 10);
    demand(n <= wt_table_quantize_N_MAX, "invalid {n}");
    demand(wia_sum <= wt_table_quantize_WIA_SUM_MAX, "invalid {wia_sum}");
    
    /* Compute the sum {wfa_sum} and max {wfa_max} of absolute input weights: */
    double wfa_sum = 0;
    double wfa_corr = 0; /* For Kahan summation. */
    double wfa_max = 0;
    for (uint32_t k = 0;  k < n; k++)
      { double wfak = fabs(wf[k]);
        demand(isfinite(wfak), "weights cannot be infinite or NAN");
        /* Kahan's summation: */
        double tmp_corr = wfak - wfa_corr;
        double tmp_sum = wfa_sum + tmp_corr;
        wfa_corr = (tmp_sum - wfa_sum) - tmp_corr;
        wfa_sum = tmp_sum;
        if (wfak > wfa_max) { wfa_max = wfak; }
      }
    if (debug) { fprintf(stderr, "    wfa_sum = %24.16e  wfa_max = %24.16e\n", wfa_sum, wfa_max); }
      
    if (wfa_sum == 0)
      { /* All weights are zero: */
        for (uint32_t k = 0;  k < n; k++) { wi[k] = 0; }
        return 0;
      }
      
    /* !!! Consider doing binsearch for {log(scale)} instead of {scale} !!! */
  
    /* Determine the max value for {scale}: */
    assert(wfa_max > 0);
    assert(wfa_sum > 0);
    double shi_elm = ((double)wt_table_quantize_WIA_MAX)/wfa_max;
    if (! isfinite(shi_elm)) { /* Gosh! Hack: */ shi_elm = 1.0e100; }
    
    double shi_sum = ((double)wt_table_quantize_WIA_SUM_MAX)/wfa_sum;
    if (! isfinite(shi_sum)) { /* Gosh! Hack: */ shi_sum = 1.0e100; }
    
    double shi = fmin(shi_elm, shi_sum);
    if (debug) { fprintf(stderr, "    shi_elm = %24.16e  shi_sum = %24.16e  shi = %24.16e\n", shi_elm, shi_sum, shi); }
    
    /* Determine the min value for {scale}, assuming all weights get rounded up: */
    double slo = ((double)wia_sum)/(wfa_sum + (double)n);
    if (debug) { fprintf(stderr, "    slo = %24.16e", slo); }
    if (slo > shi)
      { slo = shi;
        if (debug) { fprintf(stderr, " reduced to %24.16e", slo); }
      }
    if (debug) { fprintf(stderr, "\n"); }
      
    auto uint64_t quantize(double sc);
      /* Performs the scaling and quantization algorithm 
        with scale factor {sc}, which should be between
        {slo} and {shi}. */
        
    auto void debug_search
      ( uint32_t xit, 
        double xslo, uint64_t xFlo, 
        double xsmd, uint64_t xFmd,
        double xshi, uint64_t xFhi
      );
    
    /* Binary search for the best scale: */
    uint64_t Flo = quantize(slo);
    uint64_t Fhi = quantize(shi);
    if (debug) { fprintf(stderr, "    Flo = %24lu  Fhi = %24lu\n", Flo, Fhi); }
    assert(Flo <= Fhi);
    
    #define MAX_ITER 100

    double smd = NAN; /* Teh current/final guess. */
    uint64_t Fmd = 0; /* The sum of abs rounded weights for {smd}. */
    uint32_t iter = 0;
    while (TRUE)
      { if (debug) { debug_search(iter, slo, Flo, smd, Fmd, shi, Fhi); }
        
        /* At this point the final scale must be between {slo} and {shi}. */
        if (wia_sum <= Flo) { smd = slo; Fmd = quantize(smd); break; }
        if (wia_sum >= Fhi) { smd = shi; Fmd = quantize(smd); break; }
        
        /* At this point the right scale is strictly be between {slo} and {shi}. */
        if (isnan(smd)) 
          { /* Secant step: */
            double r = ((double)(wia_sum - Flo))/((double)(Fhi - Flo));
            smd = (1-r)*slo + r*shi;
            double ds = shi - slo;
            smd = fmax(slo + ds/16, fmin(shi - ds/16, smd));
          }
        else if (wia_sum == Fmd) 
          { break; }
        else
          { /* Flip-and-double step: */
            double s0 = (Fmd < wia_sum ? slo : shi); /* Move away from {s0}. */
            double s1 = (Fmd < wia_sum ? shi : slo); /* Move towards {s1}. */
            double ds = smd - s0;
            if (fabs(s1 - s0) > 4*fabs(ds))
              { smd = smd + 2*ds; }
            else
              { smd = (smd + s1)/2; }
          }
        Fmd = quantize(smd);
        iter++;
        if (wia_sum == Fmd) { break; }
        if (iter >= MAX_ITER) {  break; }
      }
    assert(! isnan(smd));
    if (debug) 
      { int64_t err = ((int64_t)Fmd) - ((int64_t)wia_sum);
        fprintf(stderr, "  final sum = %lu  error = %ld\n", Fmd, err); 
      }
    return Fmd;
    
    uint64_t quantize(double sc)
      { /* Round entries, compute actual sum of abs values: */
        uint64_t F = 0;
        for (uint32_t k = 0;  k < n; k++) 
          { double wfk = wf[k];
            if (wfk == 0) 
              { wi[k] = 0; }
            else
              { uint64_t wiak = (uint64_t)(floor(sc * fabs(wfk) + 0.5)); 
                if (keep_nz) { wiak = 1; }
                assert(wiak <= wt_table_quantize_WIA_MAX);
                wi[k] = (int32_t)(wfk < 0 ? -wiak : +wiak); 
                assert(F <= wt_table_quantize_WIA_SUM_MAX - wiak); 
                F += wiak;
              }
          }
        return F;
      }

    void debug_search
      ( uint32_t xit, 
        double xslo, uint64_t xFlo, 
        double xsmd, uint64_t xFmd,
        double xshi, uint64_t xFhi
      )
      { fprintf(stderr, "  iter = %d", xit);
        fprintf(stderr, "  min = %24.16e : %20lu", xslo, xFlo);  
        if (! isnan(smd)) 
          { fprintf(stderr, "  try = %24.16e : %20lu", xsmd, xFmd); }
        else
          { fprintf(stderr, "        %24s   %20s", "", ""); }
        fprintf(stderr, "  max = %24.16e : %20lu", xshi, xFhi); 
        fprintf(stderr, "\n");
      }
  }

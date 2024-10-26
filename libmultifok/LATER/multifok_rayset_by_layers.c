/* See {multifok_rayset.h}. */
/* Last edited on 2024-10-15 10:45:36 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

// #include <bool.h>
// #include <affirm.h>
// #include <interval.h>
// #include <r2.h>
// #include <r3.h>
// #include <rn.h>
// #include <frgb.h>
// #include <frgb_ops.h>
// #include <jsrandom.h>
// #include <jsqroots.h>
// #include <wt_table.h>
// #include <wt_table_hann.h>

#include <multifok_rayset.h>
    
  /*  The rays are organized in {NL} layers. Layer 0 has a single vertical
    ray, with {dtilt[0] = (0,0)}. Layer {kl >= 1} has {6*2^{kl-1}} rays,
    with the offsets {dtilt[kr]} equally spaced along a circle with
    center {(0,0)}.  Thus {NR} can be 1, 7, 19, 43, ... */
    
void multifok_rayset_choose_tilts_and_weights
  ( int32_t HR, 
    double zDep, 
    int32_t *NR_P,
    r2_t **dtilt_P, 
    double **wtilt_P
  )
  {
    bool_t verbose = TRUE;
    
    demand(HR >= 0, "invalid {HR}");
    demand(zDep > 0.0, "invalid {zDep}");
    
    /* Determine the number  and radii of ray layers (including the single central ray): */
    int32_t NR; /* Total number of rays, including central ray. */
    int32_t NL; /* Number of layers, including central ray. */
    double *r_layer = NULL;  /* For {kl} in {0..NL-1}, {r_layer[kl]} is the radius of layer {kl}. */
    double *w_layer = NULL;  /* For {kl} in {0..NL-1], {w_layer[kl]} is the weight of each ray in layer {kl} */
    int32_t *n_layer = NULL; /* For {kl} in {0..NL-1], {n_layer[kl]} is the number of rays in layer {kl} */
    if (isfinite(zDep)) 
      { /* Determine the number {NL} of layers beyond the central ray, and their radii: */
        /* For {kl>=1}, assumes that there are {6*2^{kl-1}} rays in circle {kl}. */
        double sigma = 2.0/zDep; /* Radius of ideal sampling distribution. */
        multifok_rayset_choose_tilt_layers(sigma, HR, &NR, &NL, &r_layer, &w_layer, &n_layer);
        assert(NL >= 2);
        assert(NR >= 2);
      }
    else
      { demand(HR == 0, "{HR} should be 0 if {zDep} is infinite");
        NL = 1;
        NR = 1;
      }
    assert(NR >= 1);
    /* Allocate the arrays: */
    r2_t *dtilt = talloc(NR, r2_t);      /* Ray directions. */
    double *wtilt = talloc(NR, double);  /* Ray weights. */
    /* Fill the central ray: */
    dtilt[0] = (r2_t){{ 0.0, 0.0 }};
    wtilt[0] = 1.0;
    if (NL >= 2)
      { assert(r_layer != NULL);
        assert(w_layer != NULL);
        /* Generate the rays by layers: */
        int32_t kr = 1; /* Number of rays generated. */
        for (int32_t kl = 1; kl < NL; kl++)
          { double rk = r_layer[kl];
            double wk = w_layer[kl];
            int32_t nk = n_layer[kl];
            double dang = 2*M_PI/nk; /* Angular step in each layer. */
            double ang = dang/2;
            for (int32_t ks = 0; ks < nk; ks++)
              { dtilt[kr].c[0] = rad*cos(ang);
                dtilt[kr].c[1] = rad*sin(ang);
                wtilt[kr] = wk;
                kr++;
                ang += dang;
              }
          }
        assert(kr == NR);
      }
      
    if (verbose)
      { fprintf(stderr, "using %d aperture rays (central plus %d layers)\n", NR, NL-1); 
        for (int32_t kr = 0; kr < NR; kr++)
          { fprintf(stderr, "  dtilt[%3d] = ( %+12.6f %+12.6f )", kr, dtilt[kr].c[0], dtilt[kr].c[1]); 
            fprintf(stderr, "  wtilt = %12.8f\n", wtilt[kr]);
          }
      }
    /* Return the results: */
    (*NR_P) = NR;
    (*dtilt_P) = dtilt;
    (*wtilt_P) = wtilt;
  }

/* Compute the weight for all rays in this layer: */
            double cr = 0.5*(1 + cos(M_PI*rad/rMax)); /* Apodizing weight. */
            cr = 1 - (1 - cr)*(1 - cr);
            double wrl = 2.0*cr/3.0;
            /* Generate {6*kl} rays on layer {kl} */
            
void multifok_rayset_choose_tilt_layers
  ( double sigma,
    int32_t HR,
    int32_t *NR_P,
    int32_t *NL_P,
    double **r_layer_P,
    double **w_layer_P,
    int32_t **n_layer_P
  )
  { 
    /* 
      The number of rays {n_layer[kl]} is 1 if {kl} is 0, otherwise it is {6*(2^(kl-1))}.
      Thus the total number of rays {NR} is 1 if {NL} is 1,
      and {1 + 6*(2^(NL-1)-1)}. The number of layers {NL} is the smallest
      value that gives {NR >= (HR+1)*(2*HR+1)}. */
      
    int32_t NL = 1;
    int32_r NR = 1;
    int32_t NR_min = (HR+1)*(2*HR+1);
    while (NR < NR_min) { NR += 6*(1 << NL); NL++; }
    fprintf(stderr, "  choose %d layers and %d rays per sample point\n", NL, NR);
    
    /* Let {F(x,y)} be a two-dimensional Gaussian PDF that is, {F(x,y) =
      A*exp(-(x^2+y^2)/4)} for some normalization constant {A}. For each
      {ku} in {0..NU}, we set {S[ku]} to an estimate of the integral of {F(x,y)}
      outside circle of radius {u(ku) = RM*ku/NU}: */
    int23_t NU = 500;
    int32_t RM = 11.0;
    double S[NU+1]; 
    /* Ignore all cponstant factors like {A,2*PI,du} for now: */
    S[NU] = 0.0;
    for (int32_t ku = NU-1; ku >= 0; ku--)
      { /* Compute midpoint {um} of ring between {u(ku)} and {u(ku+1)}: */
        double um = RM*((double)ku + 0.5)/NU;
        /* Evaluate the Gaussian at distance {um}: */
        double Fum = exp(-(um*um)/4);
        /* Multiply {Fum} by area of ring: */
        double Dk = Fum*um;
        /* Accumulate: */
        S[ku] = S[ku+1] + Dk;
      }
    /* Now normalize {S[0..NU]} so that {S[0] = 1}: */
    double Snorm = S[0];
    S[0] = 1;
    for (int32_t ku = 1; ku <= NU; ku++) { S[ku] /= Snorm; }
    /* We must have captured most of the integral: */
    assert(S[NU-1] < 1.0e-10);
    
    /* Let the /domain/ of a layer {kl} be the region of the plane that
      is assumed to be primarily sampled by the rays layer {kl}, as
      opposed to those of other layers.
      
      We will assume that the domain or layer {kl} is the region that
      closer to the circle of that layer than to other layer circles.
      That is, for {kl} in {1..NL-1}, the boundary between the domain of
      layer {kl} and that of layer {kl-1} is assumed to be the circle
      with radius {b[kl] = (r_layer[kl-1] + r_layer[kl])/2}.
      
      For completenes we define {b[0] = 0} and {b[NL] = +INF}. Then, for
      all {kl} in {0..N-1} the domain of layer {kl} is the ring between
      circles with radii {bl[kl]} and {b[kl+1]}.
      
      Thus we try to determine {bl[1..NL-1]} so that the integral of the
      Gaussian in that ring is proportional to the number of rays in that 
      ring. */
      
      double b[NL+1];
      b[NL] = +INF;
      b[0] = 0;
      double Sb = 0; 
      int32_t kub = NU;
      for (int32_t kl = NL-1; kl > 0; kl--)
        { /* At this point {Sb} is the integral in the domains of layers {kl+1..NL-1}. */
          /* And {kub} is the largest index such that {S[kub] >= Sb}. */
          /* Compute the desired integral {Sbk} in the domain of layer {kl}: */
          double Sbk = ((double)n_layer[kl])/NR;
          /* Include the domain of layer {kl} in {Sb}: */
          Sb += Sbk;
          /* Update {kub}: */
          while ((kub > 0) && (S[kub] < Sb)) { kub--; }
          /* Define {b[kl]} as the {u(kub)}: */
          b[kl] = RM*kub/NU;
        }
        
      /* Now compute the layer radii {r_layer[0..NL-1]} from {b[0..NL]}: */
      r_layer[0] = 0.0;
      for (int32_t kl = 0; kl < NL; kl++)
        { r_layer[kl] = 
        
  }

/* See {tosl_pick_planes.h} */
/* Last edited on 2024-10-06 10:04:11 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include <haf.h>

#include <tosl.h>

#include <tosl_pick_planes.h>

/* 
  SPECIALIZED PLANE GENERATORS

  The following procedures are versions of {tosl_pick_planes}
  specialized for each {type}. Each procedure fills the pre-allocated
  array {Zplane[0..NP-1]} with strictly increasing values ranging from
  {Z0} to {Z1} (inclusive both), They require {NP} to be at least 2, and
  {Z1 - Z0} to be at least {NP-1}.
    
  Unlike {tosl_pick_planes}, these procedures do not enforce any
  parity constraint on the {Z}-coordinates chosen. To get all odd
  coordinates, these procedures should be called with {Z0,Z1} scaled by
  half, and the resulting values must be converted to odd integers by
  {Zplane[ip] = 2*Zplane[ip] + 1} for all {ip}. */

void tosl_pick_planes_frac(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose);
  /* Fills {Zplane[0..NP-1]} with values as close as possible to affine interpolation between {Z0} and {Z1}. */
  
void tosl_pick_planes_sinp(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose);
  /* Fills {Zplane[0..NP-1]} with values that approximately interpolate between {Z0} and {Z1},
    but with both random and gradual deviations from uniformity. */

void tosl_pick_planes_unif(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose);
  /* Fills {Zplane[0..NP-1]} with values between {Z0} and {Z1} (inclusive) with  "fractal" or "Brownian" 
    variations. */

tosl_coord_t *tosl_pick_planes(int32_t NP, tosl_coord_t Z0, tosl_coord_t Z1, int32_t type, int32_t verbose)
  { if (verbose) { fprintf(stderr, "  > %s ----------------------------------------------------------------------\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "    NP = %d nominal range = {%+d .. %+d} type = %d\n", NP, Z0, Z1, type); }
    
    assert(NP >= 0);
    /* End values must be odd: */
    assert((Z0 & 1) == 1);
    assert((Z1 & 1) == 1);
    
    tosl_coord_t *Zplane = malloc(NP*sizeof(tosl_coord_t));
    assert(Zplane != NULL);
    
    if (NP > 0)
      { assert(Zplane != NULL);
        if (NP == 1)
          { assert(Z0 == Z1);
            Zplane[0] = Z0;
          }
        else
          { /* Cut extremes in half, rounding down even if negative: */
            tosl_coord_t Z0_half = (Z0-1)/2;
            tosl_coord_t Z1_half = (Z1-1)/2;
            switch (type)
              { case 0: tosl_pick_planes_unif(NP, Zplane, Z0_half, Z1_half, verbose); break;
                case 1: tosl_pick_planes_sinp(NP, Zplane, Z0_half, Z1_half, verbose); break;
                case 2: tosl_pick_planes_frac(NP, Zplane, Z0_half, Z1_half, verbose); break;
                default: assert(0);
              }
            /* Convert to odd integers, just doubling the spacing: */
            for (tosl_plane_id_t ip = 0; ip < NP; ip++) 
              { tosl_coord_t Zp = 2*Zplane[ip] + 1;
                /* Ensure strictly increasing, and between {Z0..Z1}: */
                if (verbose  && (NP <= 100)) { fprintf(stderr, "    Zplane[%9d] = %+9d\n", ip, Zp); }
                if (ip == 0)
                  { assert(Zp == Z0); }
                else
                  { assert(Zp > Zplane[ip-1]);
                    if (ip == NP-1)
                      { assert(Zp == Z1); }
                    else
                      { assert(Zp < Z1); }
                  }
                Zplane[ip] = Zp;
              }
          }
      }
    if (verbose) { fprintf(stderr, "  < %s ----------------------------------------------------------------------\n", __FUNCTION__); }
    return Zplane;
  }          

void tosl_pick_planes_unif(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose)
  { assert(NP >= 2);
    assert(Z1 - Z0 >= NP-1);
    for (tosl_plane_id_t ip = 0; ip < NP; ip++)
      { /* Generate a step with gradual and random variation: */
        double f = ((double)ip)/((double)NP-1); /* Relative position of plane. */
        Zplane[ip] = (tosl_coord_t)floor((1-f)*Z0 + f*Z1 + 0.5);
      }
  }

void tosl_pick_planes_sinp(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose)
  { 
    assert(NP >= 2);
    assert(Z1 - Z0 >= NP-1);

    /* The procedure computes {Zplane[ip] ≈ Z0 + A*ip + B*(sin(W*ip + D)
      - sin(D)) + C*r[ip]}, rounded to the nearest integer; where {W} is
      such that there are an integer number of {sin} cycles as {ip}
      ranges in {0..NP-1}, and {r[ip]} is arandom number in {[-1 _ +1]},
      except that {r[0] = r[NP-1] = 0}.
      
      The phase {D} can be any angle in {[0 _ 2*PI]}.  The coefficient {A}
      must be {(Z1 - Z0)/(NP-1)} so that {Zplane[NP-1]} is {Z1}.
      The coefficients {B} and {C} must be small enough so that the difference
      {dZ[ip] = Zplane[ip]-Zplane[ip-1]} is at least 1.
      
      The increment {dZ[ip]} is at least {A - B*W - 2*C}, so we must choose {B,C}
      such that {B*W + 2*C <= A - 1}.
    */
    double A = ((double)Z1-Z0)/((double)NP-1);
    assert(A >= 1.0);
    double W  = 4*M_PI; /* Freq of gradual variation. Must be an even multiple of {PI}. */
    double C = 0.05*(A - 1.0);
    double B = 0.85*(A - 1.0 - 2*C)/W; 
    double D  = ((double)random())/((double)INT32_MAX)*2*M_PI;
    if (verbose) { fprintf(stderr, "    A = %.6f B = %.6f C = %.6f\n", A, B, C); }
    
    for (tosl_plane_id_t ip = 0; ip < NP; ip++)
      { double r; /* Random perturbation term in {[-1 + +1]}. */
        if ((ip == 0) || (ip == NP-1))
          { r = 0.0; }
        else
          { r = 2*((double)random())/((double)INT32_MAX) - 1.0; }
        double Zf = Z0 + A*ip + B*(sin(W*ip + D) - sin(D)) + C*r; 
        Zplane[ip] = (tosl_coord_t)floor(Zf + 0.5);
      }
  }

void tosl_pick_planes_frac(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose)
  { 
    assert(NP >= 2);
    assert(Z1 - Z0 >= NP-1);

    auto void fill_planes(tosl_plane_id_t ipmin, tosl_plane_id_t ipmax);
      /* Fills {Zplane[ipmin+1..ipmax-1]} with values in {Zplane[ipmin]+1..Zplane[ipmax]-2}. */

    Zplane[0] = Z0;
    Zplane[NP-1] = Z1;

    fill_planes(0, NP-1);

    return;
    
    void fill_planes(tosl_plane_id_t ipmin, tosl_plane_id_t ipmax)
      { int32_t dip = ipmax - ipmin;
        assert(dip > 0);
        
        tosl_coord_t dZ = Zplane[ipmax] - Zplane[ipmin];
        assert(dZ >= dip);
        
        if (dip == 1) { /* Nothing to do: */ return; }
        tosl_plane_id_t ipctr = (ipmin + ipmax)/2; /* Middle of index range. */
        assert((ipmin < ipctr) && (ipctr < ipmax));
        if (dip == 2) { /* Not much to do: */ Zplane[ipctr] = (Zplane[ipmin] + Zplane[ipmax])/2; return; }

        /* Choose a {Zplane[ipctr]} off the affine interpolation: */
        double dZavg = ((double)dZ)/((double)dip); /* Average gap between planes. */
        
        /* Compute desired min and max gap Zmin,dZmax} between planes: */
        double fZvar = 4.0; /* Desired ratio between max and min gaps. */
        double dZmin = fmax(1.0, dZavg/sqrt(fZvar)); /* Desired min gap. */
        double dZmax = fZvar*dZmin; /* Desired max gap. */
        
        if (((random()/256) & 1) == 0)
          { /* Set {Zplane[ipctr]} as low as possible given {dZmin,dZmax}: */
            tosl_coord_t Z1min = Zplane[ipmin] + (tosl_coord_t)ceil(dZmin*(ipctr - ipmin));
            tosl_coord_t Z2min = Zplane[ipmax] - (tosl_coord_t)floor(dZmax*(ipmax - ipctr));
            Zplane[ipctr] = (Z1min > Z2min ? Z1min : Z2min);
          }
        else
          { /* Set {Zplane[ipctr]} as high as possible given {dZmin,dZmax}: */
            tosl_coord_t Z1max = Zplane[ipmin] + (tosl_coord_t)ceil(dZmax*(ipctr - ipmin));
            tosl_coord_t Z2max = Zplane[ipmax] - (tosl_coord_t)floor(dZmin*(ipmax - ipctr));
            Zplane[ipctr] = (Z1max < Z2max ? Z1max : Z2max);
          }
        
        /* Spacing should be OK: */
        assert((Zplane[ipctr] - Zplane[ipmin]) >= (ipctr - ipmin));
        assert((Zplane[ipmax] - Zplane[ipctr]) >= (ipmax - ipctr));
        
        fill_planes(ipmin, ipctr);
        fill_planes(ipctr, ipmax);
     }
  }

void tosl_pick_planes_print_stats
  ( FILE *wr,
    int32_t NP,
    tosl_coord_t Zplane[],
    tosl_coord_t Zmin,
    tosl_coord_t Zmax
  )
  {
    tosl_coord_t Z0 = Zplane[0];
    tosl_coord_t Z1 = Zplane[NP-1];

    fprintf(wr, "    selected %d planes spanning Z range = {%+d ..%+d}", NP, Z0, Z1); 
    fprintf(wr, " in {%+d ..%+d}\n", Zmin, Zmax); 
    assert(Z0 >= Zmin);
    assert(Z1 <= Zmax);
    
    /* Get max, min, avg, and dev of gaps between planes, max dev from linearity: */ 
    int32_t gap_min = INT32_MAX; /* Max gap between planes. */
    int32_t gap_max = INT32_MIN; /* Min gap between planes. */
    int32_t err_max = 0;  /* Largest deviation from linearity. */
    tosl_plane_id_t ip_err_max = -1;  /* Index of plane with largest deviation. */
    double sum_gap = 0.0;
    double sum_gap2 = 0.0; 
    tosl_coord_t Zprev = Zmin;
    for (tosl_plane_id_t ip = 0; ip <= NP; ip++)
      { tosl_coord_t Zp = (ip == NP ? Zmax : Zplane[ip]);
        int32_t gap = Zp - Zprev;
        if (ip == 0) 
          { fprintf(wr, "    gap from Zmin to Z0 = %d\n", gap); }
        else if (ip == NP) 
          { fprintf(wr, "    gap from Zplane[%d] to Zmax = %d\n", ip-1, gap); }
        else
          { if (gap < gap_min) { gap_min = gap; } 
            if (gap > gap_max) { gap_max = gap; } 
            sum_gap += (double)gap;
            sum_gap2 += gap*(double)gap;
          }
        assert(gap >= 2);
        Zprev = Zp;
        if (ip < NP)
          { tosl_coord_t Zp_lin = Z0 + (tosl_coord_t)floor(((double)ip)*(Z1 - Z0)/((double)NP-1));
            int32_t err = Zp - Zp_lin;
            if (abs(err) > abs(err_max)) { err_max = err; ip_err_max = ip; }
          }
      }
    double gap_avg = sum_gap/(NP-1);
    double gap_dev = sqrt(sum_gap2/(NP-1) - gap_avg*gap_avg);
    fprintf(wr, "    gap range = {%+d ..%+d}", gap_min, gap_max); 
    fprintf(wr, " avg = %.2f dev = %.2f\n", gap_avg, gap_dev); 
    fprintf(wr, "    largest deviation from linearity = %+d at Zplane[%d]\n", err_max, ip_err_max); 
  }

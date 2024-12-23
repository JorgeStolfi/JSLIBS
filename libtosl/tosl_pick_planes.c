/* See {tosl_pick_planes.h} */
/* Last edited on 2024-12-21 11:26:12 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include <haf.h>
#include <flt.h>
#include <jsrandom.h>

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

void tosl_pick_planes_unif(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose);
  /* Fills {Zplane[0..NP-1]} with values between {Z0} and {Z1} (inclusive) with  "fractal" or "Brownian" 
    variations. */

void tosl_pick_planes_sinp(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose);
  /* Fills {Zplane[0..NP-1]} with values that approximately interpolate between {Z0} and {Z1},
    but with both random and gradual deviations from uniformity. */

void tosl_pick_planes_frac(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose);
  /* Fills {Zplane[0..NP-1]} with values as close as possible to affine interpolation between {Z0} and {Z1}. */
  
void tosl_pick_planes_bvar(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose);
  /* Fills {Zplane[0..NP-1]} with values whose spacing varies at most by a factor of 4. */
  
tosl_coord_t *tosl_pick_planes(int32_t NP, tosl_coord_t Z0, tosl_coord_t Z1, int32_t type, int32_t verbose)
  { if (verbose) { fprintf(stderr, "  > --- %s --------------------\n", __FUNCTION__); }
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
                case 3: tosl_pick_planes_bvar(NP, Zplane, Z0_half, Z1_half, verbose); break;
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
    if (verbose) { fprintf(stderr, "  < --- %s --------------------\n", __FUNCTION__); }
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
    double W  = 4*M_PI/(NP-1); /* Freq of gradual variation. Must be an even multiple of {PI}. */
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
      /* Fills {Zplane[ipmin+1..ipmax-1]} with values in {Zplane[ipmin]+1..Zplane[ipmax]-1}. */

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

        /* Compute desired min and max gap Zmin,dZmax} between planes: */
        double dZavg = ((double)dZ)/((double)dip); /* Average gap between planes. */
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

void tosl_pick_planes_bvar(int32_t NP, tosl_coord_t Zplane[], tosl_coord_t Z0, tosl_coord_t Z1, int32_t verbose)
  { 
    bool_t debug = 0;
    
    assert(NP >= 2);
    assert(Z1 - Z0 >= NP-1);

    /* Compute desired min and max gap Zmin,dZmax} between planes: */
    tosl_coord_t dZ = Z1 - Z0;
    assert(dZ >= NP-1);
    
    double dZavg = ((double)dZ)/((double)NP-1); /* Average gap between planes. */
    double fZvar = 4.0; /* Desired ratio between max and min gaps. */
    double dZmin = fmax(1.0, dZavg/sqrt(fZvar)); /* Desired min gap. */
    double dZmax = fZvar*dZmin; /* Desired max gap. */
    if (debug) { fprintf(stderr, "  dZavg = %.4f  dZmin = %.4f  dzMax = %.4f\n", dZavg, dZmin, dZmax); }

    auto void fill_planes(tosl_plane_id_t ipmin, tosl_plane_id_t ipmax);
      /* Fills {Zplane[ipmin+1..ipmax-1]} with values in {Zplane[ipmin]+1..Zplane[ipmax]-1}. */

    Zplane[0] = Z0;
    Zplane[NP-1] = Z1;

    fill_planes(0, NP-1);

    return;
    
    void fill_planes(tosl_plane_id_t ipmin, tosl_plane_id_t ipmax)
      { if (debug) { fprintf(stderr, "  ip = %8d .. %8d  Zplane = %+12d ..%+12d\n", ipmin, ipmax, Zplane[ipmin], Zplane[ipmax]); }
        int32_t dip = ipmax - ipmin;
        assert(dip > 0);
        
        if (dip == 1)
          { /* Nothing to do: */
            assert(Zplane[ipmin] < Zplane[ipmax]);
            return;
          }
        tosl_plane_id_t ipctr = (ipmin + ipmax)/2; /* Middle of index range. */
        assert((ipmin < ipctr) && (ipctr < ipmax));
        if (dip == 2)
          { /* Not much to do: */
            assert(Zplane[ipmax] - Zplane[ipmin] >= 2);
            Zplane[ipctr] = (Zplane[ipmin] + Zplane[ipmax])/2; 
            return;
          }
        
        /* Find min and max values for {Zplane[ipctr]} given {dZmin,dZmax}: */
        assert((Zplane[ipmax] - Zplane[ipmin]) >= (ipmax - ipmin));
        
        tosl_coord_t Z1min = Zplane[ipmin] + (tosl_coord_t)floor(dZmin*(ipctr - ipmin));
        tosl_coord_t Z2min = Zplane[ipmax] - (tosl_coord_t)ceil(dZmax*(ipmax - ipctr));
        tosl_coord_t ZPmin = (Z1min > Z2min ? Z1min : Z2min);
        if (debug) { fprintf(stderr, "  Z1min = %+12d  Z2min = %+12d ZPmin = %+12d\n", Z1min, Z2min, ZPmin); }
        
        tosl_coord_t Z1max = Zplane[ipmin] + (tosl_coord_t)ceil(dZmax*(ipctr - ipmin));
        tosl_coord_t Z2max = Zplane[ipmax] - (tosl_coord_t)floor(dZmin*(ipmax - ipctr));
        tosl_coord_t ZPmax = (Z1max < Z2max ? Z1max : Z2max);
        if (debug) { fprintf(stderr, "  Z1max = %+12d  Z2max = %+12d ZPmax = %+12d\n", Z1max, Z2max, ZPmax); }
        
        /* This loop should end because {Zplane[ipmax] - Zplane[ipmin] >= ipmax - ipmin}. */
        while (ZPmin > ZPmax)
          { if (ZPmin - Zplane[ipmin] > ipctr - ipmin) { ZPmin--; }
            if (Zplane[ipmax] - ZPmax > ipmax - ipctr) { ZPmax++; }
          }
        
        /* Pick a value for {Zplane[ipctr]} in {ZPmin .. ZPmax}: */
        double f = ((double)random())/((double)INT32_MAX);
        tosl_coord_t ZPctr = ZPmin + (tosl_coord_t)floor(f*(ZPmax - ZPmin) + 0.5);
        assert((ZPmin <= ZPctr) && (ZPctr <= ZPmax));
        Zplane[ipctr] = ZPctr;
        
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
    int32_t Z_err_max = 0;  /* Largest deviation in abs of {ip -> Zp} from linearity. */
    int32_t i_err_max = 0;  /* Largest deviation in abs of {Zp -> ip} from linearity. */
    tosl_plane_id_t ip_Z_err_max = -1;  /* Index of plane with largest {ip -> Zp} deviation. */
    tosl_plane_id_t ip_i_err_max = -1;  /* Index of plane with largest {Zp -> ip} deviation. */
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
          { /* Update max linearity error.  Use ">=" instead of "=" in case it is zero. */
            tosl_coord_t Zp_lin = Z0 + (tosl_coord_t)floor(((double)ip)*(Z1 - Z0)/((double)NP-1));
            int32_t Z_err = Zp - Zp_lin;
            if (abs(Z_err) >= abs(Z_err_max)) { Z_err_max = Z_err; ip_Z_err_max = ip; }
            tosl_plane_id_t ip_lin = (tosl_plane_id_t)floor(((double)Zp-Z0)*(NP-1)/((double)Z1-Z0) + 0.5);
            int32_t i_err = ip - ip_lin;
            if (abs(i_err) >= abs(i_err_max)) { i_err_max = i_err; ip_i_err_max = ip; }
          }
      }
    double gap_avg = sum_gap/(NP-1);
    double gap_dev = sqrt(sum_gap2/(NP-1) - gap_avg*gap_avg);
    fprintf(wr, "    gap range = {%+d ..%+d}", gap_min, gap_max); 
    fprintf(wr, " avg = %.2f dev = %.2f\n", gap_avg, gap_dev); 
    fprintf(wr, "    largest deviation of {ip -> Zp} from linearity = %+d at Zplane[%d]\n", Z_err_max, ip_Z_err_max); 
    fprintf(wr, "    largest deviation of {Zp -> ip} from linearity = %+d at Zplane[%d]\n", i_err_max, ip_i_err_max); 
  }

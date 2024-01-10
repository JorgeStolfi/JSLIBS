/* See {msm_double_vec.h} */
/* Last edited on 2008-01-12 08:26:05 by stolfi */

#define msm_double_vec_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <vec.h>
#include <affirm.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <rn.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>

#include <msm_double_vec.h>

/* INTERNAL PROTOTYPES */

sign_t msm_double_vec_mutate_step(double delProb);
  /* Simulates one cycle of the natural duplication of a nucleotide 
    sequence, given its next sample {s}.
    
    With probability {delProb/2}, the procedure returns {-1},
    meaning `discard the next sample'.  
    
    With probabilty {delProb/2}, the procedure returns {+1},
    meaning `insert an extra sample before the next sample'.
    
    With the remaining probability {1-delProb}, the procedure returns {0},
    meaning `copy the next sample'. */

/* IMPLEMENTATIONS */
   
double *msm_double_vec_get_address(double_vec_t *smp, int i, bool_t circ)
  { int n = smp->ne;
    int k = (circ ? imod(i,n) : i);
    demand ((k >= 0) && (k < n), "bad index");
    return &(smp->e[k]);
  }
 
double msm_double_vec_interpolate(double_vec_t *smp, double z, bool_t circ)
  {
    double iz = floor(z);
    if (iz == z)
      { double *p = msm_double_vec_get_address(smp, (int)iz, circ);
        return (*p);
      }
    else
      { double *p0 = msm_double_vec_get_address(smp, (int)iz, circ);
        double *p1 = msm_double_vec_get_address(smp, (int)iz + 1, circ);
        double s1 = z - iz, s0 = 1 - s1;
        return s0*(*p0) + s1*(*p1);
      }
  }

void msm_double_vec_smooth(double_vec_t *smp)
  { int ns = smp->ne;
    int i;
    /* Replace the samples by their sum: */
    double s = 0;
    for (i = 0; i < ns; i++) { s += smp->e[i]; smp->e[i] = s; }
    /* Add a linear ramp to make the sequence periodic: */
    double ds = s/ns;
    for (i = 0; i < ns; i++) { smp->e[i] -= ds*i; }
  }
  
void msm_double_vec_normalize(double_vec_t *smp)
  { int ns = smp->ne;
    int i;
    /* Compute average {avg}: */
    double sum = 0;
    for (i = 0; i < ns; i++) { sum += smp->e[i]; }
    double avg = sum/ns;
    /* Shift to zero mean: */
    for (i = 0; i < ns; i++) { smp->e[i] -= avg; }
    /* Compute deviation {dev}: */
    double sum2 = 0;
    for (i = 0; i < ns; i++) { double si = smp->e[i]; sum2 += si*si; }
    double dev = sqrt(sum2/ns);
    if (dev > 0)
      { /* Scale to unit variance: */
        for (i = 0; i < ns; i++) { smp->e[i] /= dev; }
      }
  }

double_vec_t msm_double_vec_throw_normal(int ns)  
  { double_vec_t smp = double_vec_new(ns);
    rn_throw_normal(smp.ne, smp.e);
    return smp;
  } 

void msm_double_vec_mutate
  ( double_vec_t  *smpo, 
    bool_t circ,
    double mutDev, 
    double delProb,
    double_vec_t *smpd,
    msm_rung_vec_t *gv
  )
  { int no = smpo->ne; /* Number of samples in original sequence. */
    /* Allocate vectors {*smpd,*gv} slightly larger than {smpo}: */
    int nguess = 11*no/10; /* Guessed size of result. */
    (*smpd) = double_vec_new(nguess);
    (*gv) = msm_rung_vec_new(nguess);
    /* Weights for mixing samples: */
    demand((mutDev >= 0) && (mutDev <= 1), "invalid mutDev");
    double wr = mutDev; /* Weight of noise. */
    double wo = sqrt(1 - mutDev*mutDev); /* Weight of unmutated sample. */
    int nd = 0; /* Number of samples in derived sequence. */
    int ng = 0; /* Number of rungs in {gv}. */
    int io = 0; /* Index of next sample in {smpo}. */
    double sNext = (no == 0 ? 0.0 : smpo->e[imod(io,no)]);  /* Next sample from {smpo}. */
    double sPrev = (circ ? smpo->e[no-1] : sNext); /* Previous sample from {smpo}. */
    while (TRUE)
      { sign_t op = msm_double_vec_mutate_step(delProb);
        /* Check whether we are done: */
        if (op <= 0)
          { /* We will delete or copy the next sample. Does it exist? */
            if ((! circ) && (io >= no)) { break; }
            if ((circ) && (io > no)) { break; }
          }
        if (op >= 0)
          { /* Insert a new sample or copy the next sample. */
            /* Get the unmutated sample value {sd}: */
            double sPlus = (op == 00 ? sNext : (sPrev+sNext)/2);
            /* Mix it with the requested amount of noise: */
            double sr = dgaussrand();
            sPlus = wr*sr + wo*sPlus;
            /* Append it to sequence {smpd}: */
            double_vec_expand(smpd, nd); 
            smpd->e[nd] = sPlus;
            if (op == 0)
              { /* Record the inheritance of {smpd.e[nd]} from {smpo.e[io]}: */
                msm_rung_vec_expand(gv, ng); 
                gv->e[ng] = (msm_rung_t){{ io, nd }};
                ng++;
              }
            /* One more sample in {smpd}: */
            nd++;
          }
        if (op <= 0) 
          { /* Skip the sample {smp[io]}: */ 
            io++;
            if (op == 0) { /* Copied sample, update {sPrev}: */ sPrev = sNext; }
            /* Update {sNext}: */
            sNext = ((! circ) && (io >= no) ? sPrev : smpo->e[imod(io,no)]);
          }
      }
    /* Trim and return: */
    double_vec_trim(smpd, nd);
    msm_rung_vec_trim(gv, ng);
  }
    
sign_t msm_double_vec_mutate_step(double delProb)
  { 
    /* Throw a [0_1] die: */
    double toss = drandom();
    /* Note that the order of the tests is important. */
    if ((toss -= delProb/2) < 0)
      { /* Insert a new sample before the next sample: */ return +1; }
    else if ((toss -= delProb/2) < 0)
      { /* Delete the next sample: */ return -1; }
    else 
      { /* Copy the next sample with noise: */ return 0; }
  }

void msm_double_vec_write(FILE *wr, double_vec_t *smp)
  { int i;
    for (i = 0; i < smp->ne; i++)
      { fprintf(wr, "%+8.5f\n", smp->e[i]); }
    fflush(wr);
  }

void msm_double_vec_write_named(double_vec_t *smp, char *name, char *tag)
  { FILE *wr = msm_open_write(name, tag, ".seq", TRUE);
    msm_double_vec_write(wr,smp);
    fclose(wr);
  }

double_vec_t msm_double_vec_read(FILE *rd)
  { double_vec_t smp = double_vec_new(1000);
    int ns = 0; /* Counts samples (= lines) read. */
    while (TRUE)
      { double s;
        int nr = fscanf(rd, "%lf", &s);
        if (nr == EOF) { break; }
        demand(nr == 1, "format error");
        double_vec_expand(&smp, ns);
        smp.e[ns] = s;
        ns++;
      }
    /* Trim and return: */
    double_vec_trim(&smp, ns);
    return smp;
  }
  
double_vec_t msm_double_vec_read_named(char *name, char *tag)
  { FILE *rd = msm_open_read(name, tag, ".bas", TRUE);
    double_vec_t smp = msm_double_vec_read(rd);
    fclose(rd);
    return smp;
  }


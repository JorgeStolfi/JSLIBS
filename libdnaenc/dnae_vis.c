/* See dnae_vis.h */
/* Last edited on 2022-10-20 11:41:47 by stolfi */

#define dnae_vis_C_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <r3.h>
#include <argparser.h>  
#include <affirm.h>
#include <jsfile.h>

#include <msm_pairing.h>

#include <dnae_seq.h>
#include <dnae_datum.h>
#include <dnae_sample.h>

#include <dnae_vis.h>

dnae_seq_t dnae_vis_seq_read(char *name, int8_t resample_exp)
  { 
    assert(dnae_CHANNELS == 3);
    assert(resample_exp <= 0);
    
    /* Read the sequence: */
    FILE* rd = open_read(name,TRUE);
    dnae_seq_t seq = dnae_seq_read(rd);
    fclose(rd);

    if (resample_exp != 0) 
      { /* Interpolate: */
        dnae_seq_t seq_int = dnae_seq_interpolate(&seq, resample_exp); 
        dnae_seq_free_datums(&seq);
        seq = seq_int;
      }
      
    return seq;
  }

r3_vec_t dnae_vis_datums_to_points(dnae_seq_t *seq, double magnify, r3_t *pvec)
  {
    assert(dnae_CHANNELS == 3);
    int nsmp = dnae_seq_num_datums(seq);
    r3_vec_t dp = r3_vec_new(nsmp);
    int k;
    for(k = 0; k < nsmp; k++)
      { dnae_datum_t *dtk = dnae_seq_get_datum_address(seq, k);
        r3_t *dpk = &(dp.e[k]);
        int c;
        for (c = 0; c < dnae_CHANNELS; c++)
          { dpk->c[c] = dnae_sample_decode(dtk->c[c], seq->sfac.f[c]); }
        r3_scale(magnify, dpk, dpk);
        r3_add(dpk, pvec, dpk);
      }
    return dp;
  }
    
void dnae_vis_determine_visible_segments
  ( int ns, 
    dnae_seq_t seq[], 
    bool_t showMatch, 
    bool_t hideMatch, 
    double maxDist, 
    int ini[], 
    int fin[]
  )
  {
    if (maxDist >= 0)
      { if (showMatch)
          { /* Find longest internal matching seg in each sequence: */
            int minSize = 1;
            dnae_vis_find_mid_match(ns, seq, maxDist, minSize, ini, fin);
          }
        else if (hideMatch)
          { /* Find longest matching prefix and suffix: */
            int npref = dnae_vis_find_pref_suff_match(ns, seq, maxDist, +1);
            int nsuff = dnae_vis_find_pref_suff_match(ns, seq, maxDist, -1); 
            /* Show each sequence minus prefix and suffix: */
            int i;
            for (i = 0; i < ns; i++)
              { int ni = seq[i].sd.size; 
                assert((npref >= 0) && (npref <= ni));
                assert((nsuff >= 0) && (nsuff <= ni));
                if (npref + nsuff > ni + 1)
                  { /* Overlap by more than 1 datum, all effaced: */
                    ini[i] = ni;
                    fin[i] = -1;
                  }
                else
                  { /* Overlap by at most 1 datum, efface prefix and suffix minus boundary datums: */
                    ini[i] = (npref == 0 ? 0 : npref - 1);
                    fin[i] = ni - 1 - (nsuff == 0 ? 0 : nsuff - 1);
                    assert((0 <= ini[i]) && (ini[i] <= fin[i]) && (fin[i] < ni));
                  }
              }
          }
      }
    else
      { /* Show entire sequences: */
        int i;
        for (i = 0; i < ns; i++) 
          { int ni = seq[i].sd.size; 
            ini[i] = 0; 
            fin[i] = ni - 1;
          }
      }

    /* Report start and end of special section: */
    int i;
    for (i = 0; i < ns; i++) 
      { int ni = seq[i].sd.size; 
        int mi = (fin[i] >= ini[i] ? fin[i] - ini[i] + 1 : 0);
        fprintf(stderr, "sequence %d has %d effaced and %d fully-visible datums", i, ni-mi, mi);
        if (mi > 0) { fprintf(stderr, " [%d .. %d]", ini[i], fin[i]); }
        fprintf(stderr, "\n");
      } 
  }

int dnae_vis_find_pref_suff_match(int ns, dnae_seq_t seq[], double maxDist, int dir)
  { 
    bool_t debug = FALSE;
    demand((dir == -1) || (dir == +1), "invalid dir");
    double maxDist2 = maxDist*maxDist;
    /* Find the length {nmin} of the sortest sequence: */
    int nmin = (1 << 30);
    int i;
    for (i = 0; i < ns; i++) { int ni = seq[i].sd.size; if (ni < nmin) { nmin = ni; } } 
    /* Find the length {nd} of the longest matching prefix/suffix: */
    int nd = 0;
    while (nd < nmin)
      { double maxd2 = dnae_vis_max_datum_euc_distsq(ns, seq, nd, dir);
        if (debug) { fprintf(stderr, "d = %.6f  m = %.6f\n", sqrt(maxd2), maxDist); }
        if (maxd2 > maxDist2) { break; }
        nd++;
      }

    /* Report result: */
    char *which = (dir > 0 ? "prefix" : "suffix");
    fprintf(stderr, "common %s has %d datums\n", which, nd);
    return nd;
  }

void dnae_vis_find_mid_match(int ns, dnae_seq_t seq[], double maxDist, int minSize, int ini[], int fin[])
  {
    if (ns == 0) { /* Nothing to do: */ return; }
    
    double maxDist2 = maxDist*maxDist;
    int i, j;
    /* Initialize the matching segment in each sequence to be empty: */
    for (i = 0; i < ns; i++) { int ni = seq[i].sd.size; ini[i] = ni - 1; fin[i] = 0; } 

    /* Find the largest common sement in each pair, save the longest for each sequence: */
    for (i = 0; i < ns; i++)
      { int ni = seq[i].sd.size; 
        dnae_datum_scale_t *sfi = &(seq[i].sfac); 
        
        for (j = i+1; j < ns; j++)
          { /* Find the maximum matching segments between sequences {i} and {j}: */
            int nj = seq[j].sd.size; 
            dnae_datum_scale_t *sfj = &(seq[j].sfac); 
            
            /* Try all possible alignments of {seq[i]} and {seq[j]}, record the longest match: */
            
            int iniBest_i = ni;  /* Start in sequence {i} of longest {i,j} matching segment found. */
            int iniBest_j = nj;  /* Start in sequence {j} of longest {i,j} matching segment found. */
            int nmBest_ij = -1;     /* Number of datums in that matching segment: */
            
            auto void check_alignment(int ki, int kj, int64_t ng);
              /* Finds and the maximum matching subsegment in sequence {i} and {j}
                 assuming that datum {ki+r} is paired with datum {kj+r} for {r} in {0..ng-1}. */
                 
            void check_alignment(int ki, int kj, int64_t ng)
              { assert((ki >= 0) && (ki + ng <= ni));
                assert((kj >= 0) && (kj + ng <= nj));
                
                if (ng >= minSize)
                  { int rIni = 0;
                    int r;
                    for (r = 0; r < ng; r++)
                      { /* Datums {ki+rIni..ki+r-1} of {seq[i]} match {kj+rIni..kj+r-1} of {seq[j]}. */
                        dnae_datum_t *di = dnae_seq_get_datum_address(&(seq[i]), ki+r);
                        dnae_datum_t *dj = dnae_seq_get_datum_address(&(seq[j]), kj+r);
                        double d2 = dnae_datum_diffsq(di, sfi, dj, sfj);
                        if (d2 > maxDist2)
                          { /* End of run: */
                            rIni = r+1;
                          }
                        else
                          { /* Datums {rIni..r} of this alignment define a matching segment: */
                            int nm = r + rIni - 1; /* Number of datums matched. */
                            if ((nm > nmBest_ij) && (nm >= minSize)) 
                              { /* Found a better candidate for the longest matching segment: */
                                iniBest_i = ki + rIni;
                                iniBest_j = kj + rIni;
                                nmBest_ij = nm;
                              }
                          }
                      }
                  }
              }
            
            if ((ni >= minSize) && (nj >= minSize))
              { msm_pairing_enum_alignments(ni, nj, check_alignment); }
            
            /* Update the longest match for {seq[i]} and {seq[j]}: */
            if (nmBest_ij > 0)
              { if (nmBest_ij > (fin[i] - ini[i] + 1)) { ini[i] = iniBest_i; fin[i] = iniBest_i + nmBest_ij - 1; }
                if (nmBest_ij > (fin[j] - ini[j] + 1)) { ini[j] = iniBest_j; fin[j] = iniBest_j + nmBest_ij - 1; }
              }
          }
      }

    /* Report matching segments: */
    for (i = 0; i < ns; i++) 
      { int ni = (fin[i] >= ini[i] ? fin[i] - ini[i] + 1 : 0);
        if (ni == 0) 
          { fprintf(stderr, "sequence %d has no matching segment\n", i); }
        else
          { fprintf(stderr, "%d datums [%d .. %d] of sequence %d match some other sequence\n", ni, ini[i], fin[i], i);
            assert(ni >= minSize);
          }
      } 
  } 

double dnae_vis_max_datum_euc_distsq(int ns, dnae_seq_t seq[], int k, int dir)
  { 
    bool_t debug = FALSE;
    double maxd2 = -INF;
    int i, j;
    for (i = 0; i < ns; i++)
      { dnae_seq_t *seqA = &(seq[i]);
        int nsmpA = seqA->sd.size;
        int kA = (dir > 0 ? k : nsmpA - 1 - k);
        if ((kA < 0) || (kA >= nsmpA)) { return +INF; }
        dnae_datum_t *dA = dnae_seq_get_datum_address(seqA, kA);
        if (debug)
          { fprintf(stderr, "seq[%d][%d] = ", i, kA);
            dnae_datum_decoded_write(stderr, dA, &(seqA->sfac), " (", ", ", ")\n");
          }
        for (j = 0; j < i; j++) 
          { dnae_seq_t *seqB = &(seq[j]);
            int nsmpB = seqB->sd.size;
            int kB = (dir > 0 ? k : nsmpB - 1 - k);
            assert((kB >= 0) && (kB < nsmpB)); /* Since {seqB} was {seqA} before. */
            dnae_datum_t *dB = dnae_seq_get_datum_address(seqB,kB);
            double d2 = dnae_datum_euc_distsq(dA, &(seqA->sfac), dB, &(seqB->sfac));
            if (d2 > maxd2) { maxd2 = d2; }
          }
      }
    return maxd2;
  }

void dnae_vis_choose_perturbations(double perturb, int ns, r3_t pvec[])
  {
    switch(ns)
      { 
        case 0: 
          return;
        case 1:
          { pvec[0] = (r3_t){{ 0,0,0 }}; }
          return;
        case 2:
          { double d = perturb/sqrt(2.0);
            pvec[0] = (r3_t){{ -d, +d, 0.0 }};
            pvec[1] = (r3_t){{ +d, -d, 0.0 }};
          }
          return;
        case 3:
          { double d = perturb*sqrt(2.0/3.0);
            pvec[0] = (r3_t){{ -d, +d/2, +d/2 }};
            pvec[1] = (r3_t){{ +d/2, -d, +d/2 }};
            pvec[2] = (r3_t){{ +d/2, +d/2, -d }};
          }
          return;
        case 4:
          { double d = perturb/sqrt(3.0);
            pvec[0] = (r3_t){{ +d, +d, +d }};
            pvec[1] = (r3_t){{ +d, -d, -d }};
            pvec[2] = (r3_t){{ -d, +d, -d }};
            pvec[3] = (r3_t){{ -d, -d, +d }};
          }
          return;
        default:
          demand(FALSE, "not implemented");
      }
  }

void dnae_vis_parse_seqFile_options(argparser_t *pp, string_vec_t *seqFile,  int_vec_t *texture)
  {
    int nseq = 0;
    (*seqFile) = string_vec_new(10);
    (*texture) = int_vec_new(10);
    while(argparser_keyword_present(pp,"-seqFile"))
      { string_vec_expand(seqFile,nseq);
        int_vec_expand(texture,nseq);
        seqFile->e[nseq] = argparser_get_next(pp);
        if (argparser_keyword_present_next(pp,"-texture"))
          { texture->e[nseq] = (int)argparser_get_next_int(pp, 0, INT_MAX); }
        else
          { texture->e[nseq] = nseq; }
        nseq++;
      }
    string_vec_trim(seqFile,nseq);
    int_vec_trim(texture,nseq);
  }


void dnae_vis_parse_resample_option(argparser_t *pp, int8_t *resample_expP)
  {
    if (argparser_keyword_present(pp,"-resample"))
      { double st = argparser_get_next_double(pp, 1.0/((double)dna_vis_RESAMPLE_STEPS_MAX), 1.0);
        int ek = 0;
        while (st < 1) { ek--; st = st*2; }
        while (st > 1) { ek++; st = st/2; }
        if (st != 1.0) { argparser_error(pp, "initial sampling step must be a power of 2."); }
        assert(abs(ek) <= msm_seq_desc_estep_MAX);
        (*resample_expP) = (int8_t)ek;
      }
    else
      { (*resample_expP) = 0; }
  }

void dnae_vis_parse_render_options
  ( argparser_t *pp,  
    bool_t *showMatchP,
    bool_t *hideMatchP,
    double *maxDistP, 
    double *magnifyP,
    double *perturbP
  )
  {
    if (argparser_keyword_present(pp,"-hideMatch"))
      { (*hideMatchP) = TRUE; (*showMatchP) = FALSE; }
    else if (argparser_keyword_present(pp,"-hideMatch"))
      { (*hideMatchP) = FALSE; (*showMatchP) = TRUE; }
    else 
      { (*hideMatchP) = FALSE; (*showMatchP) = FALSE; }
    if ((*hideMatchP) || (*showMatchP)) 
      { (*maxDistP) = argparser_get_next_double(pp, -INF, +INF); }
    else
      { (*maxDistP) = -1.0; }

    if (argparser_keyword_present(pp,"-magnify"))
      { (*magnifyP) = argparser_get_next_double(pp, dna_vis_MAG_MIN, dna_vis_MAG_MAX); }
    else
      { (*magnifyP) = 1.0; }

    if (argparser_keyword_present(pp,"-perturb"))
      { (*perturbP) = argparser_get_next_double(pp, 0.0, 1.0); }
    else
      { (*perturbP) = dnae_vis_perturb_DEFAULT; }
  }

/* See msm_cand_vec.h */
/* Last edited on 2022-10-20 07:50:47 by stolfi */

#define msm_cand_vec_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <jswsize.h>

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>
#include <msm_cand_refine.h>
#include <msm_dyn.h>

#define msm_cand_vec_checking_level 2
  /* Define this as positive to get paranoid checking. */

#define CHKLEV(n) (msm_cand_vec_checking_level >= (n))
  /* Use {if(CHKLEV(n)){...}} for checks of paranoia level {n}. */

#define msm_cand_vec_debug_level 2
  /* Set this positive to get some debugging printouts. */

#define DEBLEV(n) (msm_cand_vec_debug_level >= (n))
  /* Use {if(DEBLEV(n)){...}} for debug printouts of level {n}. */

vec_typeimpl(msm_cand_vec_t,msm_cand_vec,msm_cand_t);

typedef double msm_subpairing_score_proc_t(int64_t r, int64_t s); 
  /* A client procedure that computes a score for a subpairing {p[r..s]}
    of some perfect pairing {p}. */

typedef void msm_subpairing_use_proc_t
  ( int64_t r,
    int64_t s,
    double score
  ); 
  /* A client procedure that uses a subpairing {p[r..s]} of some
    perfect pairing {p}, whose score is {score}. */

void msm_cand_vec_scan_alignment
  ( msm_seq_desc_t *seq0,
    int32_t i0,
    msm_seq_desc_t *seq1,
    int32_t i1,
    msm_rung_step_score_proc_t *step_score,
    int64_t minRungs,   /* Min number of rungs in a valid subpairing. */
    double minScore,    /* Ignore candidates with score lower than this. */
    msm_subpairing_use_proc_t *use  /* What to do with candidates. */
  );
  /* Scans the perfect pairing {p} between sequences {seq0,seq1} that starts
    with rungs {(i0,i1)} and has {ng} rungs. Calls {use} on some set of
    pairwise disjoint candidates found in it that have at least
    {minRungs} rungs and score at least {minScore}, hopefully with large
    scores.
    
    The score of a subpairing {p[r..s]} is assumed to be the sum of
    scores of its steps, plus the steps from nowhere to the first rung
    and from the ar rung to nowhere. */

msm_cand_vec_t msm_cand_vec_get_best_perfect
  ( msm_seq_desc_t *seq0,
    msm_seq_desc_t *seq1,
    msm_rung_step_score_proc_t *step_score,
    int64_t minRungs, 
    double minScore,
    int32_t maxCands
  )
  {
    int32_t n0 = seq0->size; 
    int32_t n1 = seq1->size; 

    msm_cand_vec_t cdv = msm_cand_vec_new(maxCands);
    int32_t ncd = 0;
    /* The current set of candidates is kept in the array
      {cdv[0..ncd-1]}. During the scan, the array is sorted
      by *decreasing* score (i.e. the worst candidate is 
      {cdv[ncd-1]}. */
      
    double lowScore = minScore; /* Current minimum score for storage. */

    auto void scan_alignment(int32_t i0, int32_t i1, int64_t ng);
      /* Scans the perfect pairing {p} that starts with rung {(i0,i1)}
        and has {ng} rungs, and adds all candidates found in it to the
        heap.
        
        A valid sub-pairing is a set of consecutive rungs from {p}
        that has at least {minRungs} and at most {min(n0,n1)} rungs. A
        candidate is a valid sub-pairing whose score is positive and
        better than that of any of proper valid sub-pairing or
        super-pairing. */

    /* Scan all maximal-length perfect pairings, collect best cands: */
    msm_pairing_enum_alignments(n0, n1,  &scan_alignment); 

    /* Trim the candidate vector to the number of cands actually found: */
    msm_cand_vec_trim(&cdv, ncd);
    
    return cdv; 

    void scan_alignment(int32_t i0, int32_t i1, int64_t ng)
      {
        auto void store_cand(int64_t kmin, int64_t kmax, double score);
          /* Stores into the sorted list {cdv} the pairing
            {(i0+k,i1+k)} for {k=kmin,..kmax}. Assumes that its
            score is {score}. The procedure is a no-op if the pairing
            is not among the best {maxCands} entries found so far, or
            is already in the list. */

        msm_cand_vec_scan_alignment(seq0, i0, seq1, i1, step_score, minRungs, lowScore, &store_cand);
     
        void store_cand(int64_t kmin, int64_t kmax, double score)
          { /* Check inex validity: */
            assert(kmin <= kmax);
            msm_rung_t gini = (msm_rung_t){{ (int32_t)(i0 + kmin), (int32_t)(i1 + kmin) }};
            msm_rung_t gfin = (msm_rung_t){{ (int32_t)(i0 + kmax), (int32_t)(i1 + kmax) }};
            assert((gini.c[0] >= 0) && (gfin.c[0] < n0));
            assert((gini.c[1] >= 0) && (gfin.c[1] < n1));
            int64_t ns = kmax - kmin + 1;
            if (DEBLEV(2))
              { fprintf
                  ( stderr, 
                    ("  got candidate from (%d,%d) to (%d,%d) steps = %" int64_d_fmt " score = %24.16e\n"), 
                    gini.c[0], gini.c[1], gfin.c[0], gfin.c[1], ns, score
                  ); 
              }
            affirm(score >= lowScore, "bad score");
            msm_pairing_t *pr = msm_pairing_perfect(gini, (int32_t)ns);
            msm_cand_t cd = msm_cand_from_pairing(seq0, seq1, pr, score);
            double frac = +INF;
            msm_cand_vec_insert(&cdv, &ncd, (int32_t)minRungs, maxCands, &cd, frac);
            /* Update the current min useful score: */
            lowScore = (ncd < maxCands ? minScore : cdv.e[ncd-1].score);
            if (DEBLEV(1)) { fprintf(stderr, "  %d candidates, min score = %24.16e\n", ncd, lowScore); }
          }
      }
  }

void msm_cand_vec_scan_alignment
  ( msm_seq_desc_t *seq0,
    int32_t i0,
    msm_seq_desc_t *seq1,
    int32_t i1,
    msm_rung_step_score_proc_t *step_score,
    int64_t minRungs,   /* Min number of rungs in a valid subpairing. */
    double minScore,    /* Ignore candidates with score lower than this. */
    msm_subpairing_use_proc_t *use  /* What to do with candidates. */
  )
  {
    int32_t n0 = seq0->size; 
    int32_t n1 = seq1->size;

    affirm((i0 >= 0) && (i0 < n0), "bad {i0}");
    affirm((i1 >= 0) && (i1 < n1), "bad {i1}");
    affirm((i0 == 0) || (i1 == 0), "bad {i0,i1}");
    
    /* Get max number of rungs {ng} with this alignment: */
    int64_t ng = imin(n0 - i0, n1 - i1);  
    if (DEBLEV(1))
      { fprintf
          ( stderr, ("\nalignment (%d,%d:%" int64_d_fmt ") minScore = %12.6f\n"),
            i0, i1, ng, minScore
          );
      }
    affirm(ng >= 0, "bad {ng}");

    /* Trivial cases: */
    if (ng < minRungs) { return; }
    
    /* Let {p} be the perfect pairing that starts with rung {i0,i1} and
      has {ng} rungs. Uses dynamic programming in one dimension to find
      the sub-pairing of {p} with maximum score that ends with each rung
      {k}. */

    double psc[ng]; /* {psc[k]} is the score of the best sub-pairing that ends at rung {k}. */
    
    /*  Scans the rungs {(i0+k,i1+k)} for {k} from {0} to {ng-1}. At
      each rung there are two choices: start a new sub-pairing with that
      rung, or extend the current one to that rung. In the former case
      the current rung is terminated and reported through {use}. */
    
    auto void use_previous(int64_t kini, int64_t kmax);
      /* Report the best candidate that starts with rung {kini}, 
         ends at or before rung {kmax},and has at least {minRungs} rungs,
         provided its score is at least {minScore}. */
    
    int64_t kcur = 0; /* First rung of current pairing. */
    msm_rung_t h = msm_rung_none;  /* Previous rung. */

    int64_t k;
    for (k = 0; k < ng; k++)
      { 
        msm_rung_t g = (msm_rung_t){{ (int32_t)(i0 + k), (int32_t)(i1 + k) }};
        
        double score_new = step_score(seq0, seq1, NULL, &g); /* Score for starting a new sub-pairing at rung {k}. */
        double score_old = ( k <= 0 ? -INF : psc[k-1] + step_score(seq0, seq1, &h, &g) ); /* Score for extending. */
        
        if (score_new > score_old)
          { /* Better start a new candidate here: */
            use_previous(kcur, k - 1);
            psc[k] = score_new;
            /* Remember candidate start: */
            kcur = k;
          }
        else
          { /* Extend old candidate: */
            psc[k] = score_old;
          }
      }
    use_previous(kcur, ng-1);
    return;
    
    void use_previous(int64_t kini, int64_t kmax)
      { 
        /* Find the best sub-pairing that ends at or before rung {kmax}: */
        int64_t kmin = kini + minRungs - 1;
        if (kmax < kmin) { /* Too short: */ return; }
        int64_t kopt = kmax;
        int64_t k = kmax;
        while (k > kmin) 
          { msm_rung_t h = (msm_rung_t){{ (int32_t)(i0 + k), (int32_t)(i1 + k) }};
            psc[k] += step_score(seq0, seq1, &h, NULL);
            if ((k < kopt) && (psc[k] > psc[kopt])) { kopt = k; }
            k--; 
          }
        /* Get the best rung: */
        double score = psc[kopt];
        if (score >= minScore) { use(kini, kopt, score); }
      }
  }

void msm_cand_vec_throw
  ( int32_t ntr,
    msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    int32_t minlen,
    int32_t maxlen,
    double atomProb,
    double diagProb,
    double skipProb,
    msm_cand_vec_t *cdv,
    int32_t *ncdP
  )
  { int32_t ncd = (*ncdP);
    int32_t i; 
    for (i = 0; i < ntr; i++) 
      { /* Choose between atomic or not: */
        bool_t atomP = (drandom() < atomProb);
        /* Choose between near-diagonal and random: */
        bool_t diagP = (drandom() < diagProb);
        /* Choose length of candidate: */
        int32_t len = int32_abrandom(minlen, maxlen);
        /* Generate random candidate: */
        msm_cand_vec_expand(cdv, ncd);
        cdv->e[ncd] = msm_cand_throw(seq0, seq1, len, atomP, diagP, 0.10);
        ncd++;
      }
    (*ncdP) = ncd;
  }

void msm_cand_vec_free(msm_cand_vec_t *cdv)
  {
    int32_t i;
    for (i = 0; i < cdv->ne; i++) { msm_cand_free(&(cdv->e[i])); }
    free(cdv->e);
  }

#define msm_cand_vec_type_name "msm_cand_vec"
#define msm_cand_vec_old_format "2008-01-11"
#define msm_cand_vec_format "2013-10-16"

void msm_cand_vec_write(FILE *wr, msm_cand_vec_t *cdv)
  { int32_t ncd = cdv->ne;
    /* Find maximum abs of {id,strlen(name),index}: */
    int32_t idmax = 1; 
    int32_t nlenmax = 1;
    int32_t ixmax = 1;
    int32_t i;
    for (i = 0; i < ncd; i++)
      { msm_cand_t *cd = &(cdv->e[i]);
        int32_t j;
        for (j = 0; j < 2; j++)
          { msm_seq_desc_t *sj = &(cd->seq[j]);
            int32_t idj = (sj->id >= 0 ? sj->id : 10*(-sj->id)); 
            if (idj > idmax) { idmax = idj; }
            int32_t lenj = (sj->name == NULL ? 0 : (int32_t)strlen(sj->name));
            if (lenj > nlenmax) { nlenmax = lenj; }
            int32_t ixj = (sj->size <= 0 ? 0 : sj->size - 1);
            if (ixj > ixmax) { ixmax = ixj; }
          }
      }
    /* Compute the max number of digits needed for {i0,name,index}: */
    int32_t idSize = digits(idmax);
    int32_t nameSize = nlenmax;
    int32_t ixSize = digits(ixmax);
    filefmt_write_header(wr, msm_cand_vec_type_name, msm_cand_vec_format);
    fprintf(wr, "ncands = %d\n", ncd);
    for (i = 0; i < ncd; i++)
      { msm_cand_write(wr, NULL, &(cdv->e[i]), idSize, nameSize, ixSize, NULL,TRUE);
        fputc('\n', wr);
      }
    filefmt_write_footer(wr, msm_cand_vec_type_name);
    fflush(wr);
  }

void msm_cand_vec_write_named(msm_cand_vec_t *cdv, char *name, char *tag, char *ext)
  { FILE *wr = msm_open_write(name, tag, ext, TRUE);
    msm_cand_vec_write(wr, cdv);
    fclose(wr);
  }
  
msm_cand_vec_t msm_cand_vec_read(FILE *rd)
  { /* Check and skip the file header: */
    filefmt_read_header(rd, msm_cand_vec_type_name, msm_cand_vec_format);
    /* Skip comment lines, if any: */
    (void)filefmt_read_comment(rd, '|');
    /* Read the header fields: */
    bool_t ncd = nget_int32(rd, "ncands"); fget_eol(rd);
    /* Read the candidates: */
    msm_cand_vec_t cdv = msm_cand_vec_new(ncd);
    int32_t k; 
    for (k = 0; k < ncd; k++)
      { cdv.e[k] = msm_cand_read(rd, NULL, NULL); 
        fget_eol(rd);
      }
    /* Check and skip file footer: */
    filefmt_read_footer(rd, msm_cand_vec_type_name);
    return cdv;
  }
 
msm_cand_vec_t msm_cand_vec_read_named(char *name, char *tag, char *ext)
  { FILE *rd = msm_open_read(name, tag, ext, TRUE);
    msm_cand_vec_t cdv = msm_cand_vec_read(rd);
    fclose(rd);
    return cdv;
  }

msm_cand_vec_t msm_cand_vec_map(msm_cand_vec_t *cdv, msm_seq_desc_t *s0_new, msm_seq_desc_t *s1_new)
  { fprintf(stderr, "  mapping %d candidates ...\n", cdv->ne);
    msm_cand_vec_t cdvnew = msm_cand_vec_new(cdv->ne);
    int32_t ic;
    for (ic = 0; ic < cdv->ne; ic++)
      { cdvnew.e[ic] = msm_cand_map(&(cdv->e[ic]), s0_new, s1_new); }
    fprintf(stderr, "  mapped %d candidates.\n", cdvnew.ne);
    return cdvnew;
  }
 
msm_cand_vec_t msm_cand_vec_interpolate(msm_cand_vec_t *cdv)
  { msm_cand_vec_t cdvnew = msm_cand_vec_new(cdv->ne);
    int32_t ic;
    for (ic = 0; ic < cdv->ne; ic++)
      { cdvnew.e[ic] = msm_cand_interpolate(&(cdv->e[ic])); }
    return cdvnew;
  }
 
msm_cand_vec_t msm_cand_vec_make_increasing(msm_cand_vec_t *cdv, int32_t minIncrEach, int32_t minIncrSum)
  { msm_cand_vec_t cdvnew = msm_cand_vec_new(cdv->ne);
    int32_t ic;
    for (ic = 0; ic < cdv->ne; ic++)
      { cdvnew.e[ic] = msm_cand_make_increasing(&(cdv->e[ic]), minIncrEach, minIncrSum); }
    return cdvnew;
  }

msm_cand_vec_t msm_cand_vec_refine
  ( msm_cand_vec_t *cdvold,
    msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    int32_t delta,
    int32_t kappa, 
    int32_t expand,
    int32_t shrink,
    int32_t maxUnp, 
    msm_rung_step_score_proc_t *step_score,
    bool_t verbose,
    msm_dyn_tableau_t *tb, 
    int32_t minCover,
    int32_t maxCands,
    double frac
  )
  { msm_cand_vec_t cdvnew = msm_cand_vec_new(cdvold->ne);
    int32_t ic;
    int32_t ncnew = 0; /* Number of candidates that were kept. */
    int32_t n_entries = 0;
    int32_t n_steps = 0;
    fprintf(stderr,"refining candidate vector...\n");
    for (ic = 0; ic < cdvold->ne; ic++)
      { msm_cand_t *cdoi = &(cdvold->e[ic]);
        demand(msm_seq_desc_equal(&(cdoi->seq[0]), seq0, FALSE), "wrong seq 0");
        demand(msm_seq_desc_equal(&(cdoi->seq[1]), seq1, FALSE), "wrong seq 1");
        
        msm_cand_t cdni = msm_cand_refine
          ( cdoi, delta, kappa, expand, shrink, maxUnp, step_score, verbose, tb, &n_steps, &n_entries );
          
        msm_cand_vec_insert(&cdvnew, &ncnew, minCover, maxCands, &cdni, frac);
      }
    msm_cand_vec_trim(&cdvnew, ncnew);
    fprintf
      ( stderr, "refined candidate vactor: %d entries, %.1f steps per entry\n",
        n_entries, ((double)n_steps)/n_entries
      );
    return cdvnew;
  }

void msm_cand_vec_insert
  ( msm_cand_vec_t *cdv,
    int32_t *ncandP,
    int32_t minCover,
    int32_t maxCands,
    msm_cand_t *cd,
    double frac
  )
  {
    auto int32_t bubble_up(int32_t j);
      /* Assumes that the list is OK except that {cdv[j]} may be worse
        than {cdv[j-1]}. Swaps the two. */
    
    auto void delete(int32_t j);
      /* Deletes the entry {cdv[j]}, displacing all subsequent entries. 
        Decrements {ncand}. */
    
    /* Check minimum coverage requirement: */
    int32_t j;
    for (j = 0; j < 2; j++)
      { if (msm_pairing_span(cd->pr, j) < minCover) 
        { /* Candidate does not cover enough of sequence {j}: */
          msm_cand_debug("N", cd);
          return;
        }
      }
    
    int32_t k; /* Final position of inserted candidate, or {(*ncandP)} if not inserted. */

    /* Find the smallest {k} such that {cdv[k]} is definitely worse than {cd}: */
    k = (*ncandP);
    while ((k > 0) && (cdv->e[k-1].score < cd->score)) { k--; }

    if ((frac > 0.0) && (frac <= 1))
      { 
        /* Discard {cd} if there is any better or equal cand that contains or is similar to {cd}: */
        for (j = 0; j < k; j++)
          { int32_t ngnew = msm_pairing_num_rungs(cd->pr); /* Rungs in {cd}. */
            int32_t ngboth = msm_cand_count_common_rungs(&(cdv->e[j]), cd); /* Shared rungs. */
            if (ngboth >= frac*ngnew) 
              { /* Discard {cd}: */
                msm_cand_debug("C", cd);
                return;
              }
          }

        /* Remove any candidates that are worse than {cd} but contained in or similar to {cd}: */
        for (j = (*ncandP)-1; j >= k; j--)
          { int32_t ngold = msm_pairing_num_rungs(cdv->e[j].pr); /* Rungs in {cdv->e[j]}. */
            int32_t ngboth = msm_cand_count_common_rungs(&(cdv->e[j]), cd); /* Shared rungs. */
            if (ngboth >= frac*ngold) 
              { /* Discard candidate {cvd.e[j]}: */
                //msm_cand_debug("S", &(cdv->e[j]));
                delete(j);
              }
          }
      }

    if ((*ncandP) >= maxCands)
      { k = (*ncandP)-1;
        if (cd->score < cdv->e[k].score)
          { /* Discard {cd}: */
            msm_cand_debug("B", cd);
            return;
          }
        else
          { /* Delete worst candidate: */
            msm_cand_debug("-", &(cdv->e[k]));
            delete(k);
          }
      }

    /* Add {cd} to {cdv[0..ncand-1]}: */
    assert((*ncandP) < maxCands);
    k = (*ncandP);
    msm_cand_vec_expand(cdv, k);
    cdv->e[k] = (*cd);
    (*ncandP)++;
    k = bubble_up(k);
    assert((k >= 0) && (k < (*ncandP)));
    msm_cand_debug("+", &(cdv->e[k]));

    int32_t bubble_up(int32_t j)
      { while (j > 0)
          { /* Get parent: */
            int32_t k = j-1; 
            /* If {cdv[k]} is not worse than {cdv[j]}, stop: */
            if (cdv->e[k].score >= cdv->e[j].score) { break; }
            /* Swap parent with {cdv[j]}: */
            msm_cand_t pp = cdv->e[k]; cdv->e[k] = cdv->e[j]; cdv->e[j] = pp;
            j = k;
          }
        return j;
      }

    void delete(int32_t j)
      { int32_t i;
        for (i = j+1; i < (*ncandP); i++) { cdv->e[i-1] = cdv->e[i]; }
        (*ncandP)--;
      }
      
  }

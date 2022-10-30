#define PROG_NAME "msm_test_020_refining"
#define PROG_DESC "test of candidate refinement routines"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-30 11:19:08 by stolfi */

#define msm_test_020_refining_C_COPYRIGHT \
  "Copyright � 2006  by the State University of Campinas (UNICAMP)"

/* !!! Add subsampling */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -seqLength {SEQ_LENGTH} \\\n" \
  "  -expSub {EXPSUB} \\\n" \
  "  -mutDev {MUT_DEV} \\\n" \
  "  -delProb {INS_PROB} \\\n" \
  "  [ -expand {EXPAND} ] [ -shrink {SHRINK} ] \\\n" \
  "  [ -delta {DELTA} ] [ -kappa {KAPPA} ] [ -maxUnp {MAXUNP} ] \\\n" \
  "  -nCands {N_CANDS} \\\n" \
  "  [ -repeat {N_TIMES} ] \\\n" \
  "  {OUT_NAME}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program generates {N_CANDS} random candidate pairings" \
  " between two random sequences with approximately {SEQ_LENGTH}" \
  " real-valued samples, and then refines them, {N_TIMES} times.\n" \
  "\n" \
  "  The two sequences are generated by making two copies of" \
  " a random sequence, inserting samples" \
  " at random in each copy, and adding independent" \
  " Gaussian noise to each sample.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  All output files will have names starting with {OUT_NAME}.\n" \
  "\n" \
  "  All candidates are plotted together as image files" \
  " \"{OUT_NAME}-{XXX}-cd.pgm\" where {XXX} is \"ini\" (before" \
  " refinement) or \"fin\" (after it)." \
  "\n" \
  "  The sample-to-sample distance matrix for the two sequences" \
  " is plotted as file \"{OUT_NAME}-eq.pgm\"." \
  "\n" \
  "  The dynamic programming tableau for each candidate" \
  " and each iteration of the refinement procedure is plotted as image" \
  " file \"{OUTNAME}-{NNNNN}-{T}.pgm\" where {NNNNN} is a" \
  " five-digit candidate number, and {T} is a single-digit" \
  " iteration number.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -seqLength {SEQ_LENGTH}\n" \
  "    This mandatory argument specifies the" \
  " approximate length (number of samples)" \
  " of the two test sequences.\n" \
  "\n" \
  "  -expSub {EXPSUB}\n" \
  "    This mandatory argument specifies the" \
  " subsampling factor for pairing.  Rung coordinates will be" \
  " implicitly divided by {2^EXPSUB} to obtain the sample indices.\n" \
  "\n" \
  "  -mutDev {MUT_DEV}\n" \
  "    This mandatory argument specifies the" \
  " standard deviation of the Gaussian noise to be added to each " \
  " sample in each sequnce.\n" \
  "\n" \
  "  -delProb {DEL_PROB}\n" \
  "    This mandatory argument specifies the sample deletion" \
  " or insertion probability used when generating the" \
  " two sequences.\n" \
  "\n" \
  "  -delta {DELTA}\n" \
  "    This optional argument specifies the" \
  " amount of adjustment allowed for the X and Y coordinates of" \
  " internal rungs of each pairing.  The default is 3.\n" \
  "\n" \
  "  -kappa {DELTA}\n" \
  "    This optional argument specifies the" \
  " amount of X and Y extension allowed" \
  " at either end of each pairing.  The default is 6.\n" \
  "\n" \
  "  -maxUnp {MAXUNP}\n" \
  "    This optional argument specifies the" \
  " maximum unpaired datums between any two rungs of each" \
  " pairing.  The default is 6.\n" \
  "\n" \
  "  -expand {EXPAND}\n" \
  "    This optional argument specifies by how much the {R}-range" \
  " of the refined candidates may extend beyond the" \
  " original candidate's {R}-range, in both directions.  If omitted, or" \
  " if {EXPAND} is zero, the refined" \
  " {R}-range will be a subset of the orginal range.\n" \
  "\n" \
  "\n" \
  "  -shrink {SHRINK}\n" \
  "    This optional argument specifies by how much the {R}-range of" \
  " the refined candidates may shrink into the original" \
  " candidate's {R}-range, in each direction.  If omitted, or" \
  " if {SHRINK} is zero, the refined {R}-range" \
  " will be a superset of the original range.\n" \
  "\n" \
  "  -nCands {N_CANDS}\n" \
  "    This mandatory argument specifies how many candidates" \
  " to generate.\n" \
  "\n" \
  "  -repeat {N_TIMES}\n" \
  "    This optional argument specifies how many times each" \
  " candidate should be refined.  The default is 1.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  msm_test_120_mapping(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 21/dec/2006 by J. Stolfi.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " msm_test_020_refining_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <affirm.h>
#include <float_image.h>
#include <uint16_image.h>
#include <jsmath.h>
#include <argparser.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>
#include <msm_rung.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>
#include <msm_cand_refine.h>
#include <msm_dyn.h>
#include <msm_image.h>
#include <msm_image_tools.h>
#include <msm_double_vec.h>
#include <msm_test_tools.h>

typedef int32_t msm_seq_pos_t;
/* The numerator of a fractional index into in a sequence. 
   The denominator should be specified separately. */

typedef struct msm_options_t 
  { /* Sequence generation parameters: */
    int32_t seqLength;       /* Expected number of positions on each sequence. */
    double mutDev;       /* Mutation probability. */
    double delProb;      /* Insertion/deletion probability. */
    int32_t expSub;          /* Exponent of subsampling factor (denominator for {msm_seq_pos_t}s). */
    /* Candidate generation parameters: */
    int32_t nCands;          /* Number of candidates to generate. */
    /* Refinement iteration parameter: */
    int32_t repeat;          /* How many times to refine each candidate. */
    /* Refinement parameters: */
    int32_t delta;           /* Half-width of tableau around original pairing. */
    int32_t kappa;           /* Extension of tableau beyond ends of pairing. */
    int32_t expand;          /* How much the refined {R}-range may expand beyond original {R}-range. */
    int32_t shrink;          /* How much the refined {R}-range may shrink into the original {R}-range. */
    int32_t maxUnp;          /* Maximum unpaired datums allowed between consecutive rungs. */
    /* Output parameters: */
    char *outName;  /* Output file name prefix (minus extensions). */
  } msm_options_t;
  
int32_t main(int32_t argc, char**argv);

msm_options_t *msm_get_options(int32_t argc, char**argv);
  /* Parses the command line options, packs 
    them into a {msm_options_t} record. */

void msm_test_seq_make_ancestral
  ( int32_t ns, 
    msm_seq_id_t id, 
    char *name,
    int32_t expSub, 
    msm_seq_desc_t *seq, 
    double_vec_t *smp
  );
  /* Creates a test sequence with {ns} random real samples,
    subsampled by a factor {den=2^expSub}. Returns the
    sequence descriptor in {*seq} and the sample vector in {*smp}
    (which is allocated by the procedure).
    
    The sequence descriptor {*seq} will have the given {id,name} attributes,
    {estep=expSub}, {skip=0}, and {size} equal to the number of subsampling
    points, namely {(ns-1)*den+1}. Its indices will
    therefore range in {0..size-1}. */

void msm_test_seq_make_derived
  ( msm_seq_desc_t *seqo,
    double_vec_t *smpo,
    double mutDev,
    double delProb,
    msm_seq_id_t idd, 
    char *named,
    msm_seq_desc_t *seqd,
    double_vec_t *smpd,
    msm_rung_vec_t *gv
  );
  /* Creates a test sequence by mutating a given sequence. Assumes
    that the original sequence has descriptor {*seqo} and 
    samples {smpo.e[0..nso-1]}, where {nso == smpo.ne}.
    
    Returns the sequence descriptor in {*seqd}, and the samples of the
    new sequence in {*smpd} (which is allocated by the procedure).
    Its length {nsd=smpd.ne} will be similar to but not always equal to {nso}.
    
    The new sequence descriptor {*seqd} will have the given
    {id,name}, the same {rev,estep,skip} attributes as {*seqo}, and
    {size=(nso-1)*den+1-2*skip} where {den = 2^(estep)}.
    
    Also returns in {gv} a list of rungs that connect the original
    subsamples of the ancestral sequence to their copies in the new
    sequence. */

void msm_fake_initial_cands
  ( int32_t ntr, 
    msm_seq_desc_t *seq0,
    msm_seq_desc_t *seq1,
    msm_cand_vec_t *cdv,
    int32_t *ncdP
  );
  /* Generates {ntr} candidates for the two {msm_seq_desc_t}s.
    At east one candidate is near the diagonal. */
  
double msm_test_get_sample(msm_seq_desc_t *seq, double_vec_t *smp, msm_seq_pos_t i);
  /* Gets value of sample {i} from sequence {seq} whose samples are {smp}. */

double msm_test_sample_diffsq
  ( msm_seq_desc_t *seq0,
    double_vec_t *smp0,
    msm_seq_pos_t i0,
    msm_seq_desc_t *seq1,
    double_vec_t *smp1, 
    msm_seq_pos_t i1,
    bool_t debug
  );
  /* Returns the square of the difference between sample {smp0(i0/den0)}
    and {smp1(i1/den1)}. Uses linear interpolation.  */

double msm_test_step_score
  ( msm_seq_desc_t *seq0,
    double_vec_t *smp0,
    msm_seq_desc_t *seq1,
    double_vec_t *smp1,
    msm_rung_t *g, 
    msm_rung_t *h,
    double EQL,
    double DIF,
    double BRK,
    double SKP
  );
  /* Returns a numeric quality score for the step {g-->h} on sequences
    {seq0,seq1}, whose original sample vectors are {smp0,smp1}.
    
    The parameters {EQL,DIF,BRK,SKP} are the scoring points for for
    equality squared, difference squared, imperfect steps, and unpaired
    datums, respectively. Usually {EQL} is positive and the rest is
    negative. */

void msm_show_tableau
  ( char *outName,         /* File name prefix. */
    int32_t ic,                /* Candidate index. */
    int32_t ir,                /* Refinement iteration index. */
    msm_cand_t *cd,        /* Refined candidate. */
    msm_dyn_tableau_t *tb  /* Dynamic programming tableau. */
  );
  /* Writes an image called "{outName}-{NNNNN}-{R}.ppm" containing a picture
    of the tableau {tb}, with the candidate {cd} drawn on top of it. */ 

int32_t main(int32_t argc, char**argv)
  { 
    msm_options_t *o = msm_get_options(argc, argv);
    
    fprintf(stderr, "generating and writing the ancestral sequence ...\n");
    double_vec_t smpa;
    msm_seq_desc_t seqa;
    msm_seq_id_t ida = 0;
    char *namea = "A";
    msm_test_seq_make_ancestral(o->seqLength, ida, namea, o->expSub, &seqa, &smpa);
    msm_test_seq_write_and_plot_named(&seqa, &smpa, "Ancestral", o->outName, "-org", 10.0, TRUE);

    fprintf(stderr, "generating and writing the mutated copies ...\n");
    double_vec_t smp[2];
    msm_seq_desc_t seq[2];
    msm_rung_vec_t gvd[2];
    int32_t j;
    for (j = 0; j < 2; j++)
      { msm_seq_id_t seqId = j+1;
        char *seqName = (j == 0 ? "X" : "Y");
        char *title = (j == 0 ? "Derived X" : "Derived Y");
        char *tag = (j == 0 ? "-mux" : "-muy");
        msm_test_seq_make_derived
          ( &seqa, &smpa,
            o->mutDev, o->delProb, 
            seqId, seqName,
            &(seq[j]), &(smp[j]), &(gvd[j])
          );
        msm_test_seq_write_and_plot_named(&(seq[j]), &(smp[j]), title, o->outName, tag, 10.0, TRUE);
      }
    
    msm_rung_vec_t gv = msm_rung_vec_join(&(gvd[0]), 0, &(gvd[1]), 0);
    
    auto double rung_score(msm_seq_desc_t *s0, msm_seq_desc_t *s1, msm_rung_t *g);

    double rung_score(msm_seq_desc_t *s0, msm_seq_desc_t *s1, msm_rung_t *g)
      { 
        (void)msm_seq_desc_same_orig_seq(s0, &(seq[0]), TRUE);
        (void)msm_seq_desc_same_orig_seq(s1, &(seq[1]), TRUE);
        return msm_test_sample_diffsq
          ( s0, &(smp[0]), g->c[0], 
            s1, &(smp[1]), g->c[1],
            FALSE
          );
      }

    fprintf(stderr, "plotting distance matrix for test seq pair ...\n");
    msm_pairing_t *pr = msm_pairing_from_rung_vec(&gv);
    msm_image_seq_seq_score_write_named(&(seq[0]), &(seq[1]), rung_score, pr, o->outName, "-eq");
    
    fprintf(stderr, "creating initial candidates ...\n");
    msm_cand_vec_t cdvraw = msm_cand_vec_new(o->nCands);
    int32_t ncd = 0; /* Candidates are {cdvraw.e[0..ncd-1]}. */
    /* Candidate 0 is the "true" pairing: */
    msm_rung_vec_t gvtrue = msm_rung_vec_make_increasing(&gv, 1, 1);
    gvtrue = msm_rung_vec_interpolate(&gvtrue);
    msm_pairing_t *ptrue = msm_pairing_from_rung_vec(&gvtrue);
    cdvraw.e[ncd] = msm_cand_from_pairing(&(seq[0]), &(seq[1]), ptrue, 0.0);
    ncd++;
    /* Additional candidates: */
    msm_fake_initial_cands(o->nCands-1, &(seq[0]), &(seq[1]), &cdvraw, &ncd);
    msm_cand_vec_trim(&cdvraw, ncd);

    fprintf(stderr, "writing initial candidates ...\n");
    msm_cand_vec_write_named(&cdvraw, o->outName, "-ini-cd", ".cdv");
    msm_image_cand_vec_write_named(&cdvraw, &(seq[0]), &(seq[1]), FALSE, o->outName, "-ini-cd");
    
    /* Weights, suitable for one-dimensional signals: */
    double EQL = +1.0;
    double DIF = -1.0;
    double BRK = 0.00;
    double SKP = 0.00;
      
    auto double step_score(msm_seq_desc_t *s0, msm_seq_desc_t *s1, msm_rung_t *g, msm_rung_t *h);

    double step_score(msm_seq_desc_t *s0, msm_seq_desc_t *s1, msm_rung_t *g, msm_rung_t *h)
      { return msm_test_step_score(s0, &(smp[0]), s1, &(smp[1]), g, h, EQL, DIF, BRK, SKP); }

    fprintf(stderr, "refining candidates ...\n");
    msm_dyn_tableau_t tb = msm_dyn_tableau_new(); /* Dynamic programming tableau. */
    msm_cand_vec_t cdvref = msm_cand_vec_new(ncd); /* Refined cands will be {cdvref[0..ncd-1]}. */
    int32_t ic;
    for (ic = 0; ic < ncd; ic++)
      { msm_cand_t cd = cdvraw.e[ic];
        /* Refine candidate {o->repeat} times: */
        int32_t ir;
        for (ir = 0; ir < o->repeat; ir++)
          { int32_t n_steps = 0;
            int32_t n_entries = 0;
            bool_t verbose = TRUE;
            msm_cand_t cdref = msm_cand_refine
              ( &cd,
                o->delta, o->kappa, o->expand, o->shrink, o->maxUnp,
                &step_score,
                verbose,
                &tb,
                &n_steps, &n_entries
              );
            /* Display the matrix: */
            msm_show_tableau(o->outName, ic, ir, &cdref, &tb);
           /* Prepare for next iteration: */
            cd = cdref;
          }
        cdvref.e[ic] = cd;
      }
      
    fprintf(stderr, "writing final candidates ...\n");
    msm_cand_vec_write_named(&cdvref, o->outName, "-fin-cd", ".cdv");
    msm_image_cand_vec_write_named(&cdvref, &(seq[0]), &(seq[1]), FALSE, o->outName, "-fin-cd");

    return 0;
  }

void msm_test_seq_make_ancestral
  ( int32_t ns, 
    msm_seq_id_t id, 
    char *name,
    int32_t expSub, 
    msm_seq_desc_t *seq, 
    double_vec_t *smp
  )
  { /* Generate a vector of random numbers: */
    *smp = msm_double_vec_throw_normal(ns);
    /* Smooth it a few times: */
    int32_t nsmooth = 1;
    int32_t j;
    for (j = 0; j < nsmooth; j++)
      { msm_double_vec_smooth(smp); 
        msm_double_vec_normalize_avg_dev(smp);
      }
    /* Compute the subsampling factor {den}: */
    assert((expSub >= 0) && (expSub <= 10));
    int32_t den = (1 << expSub);
    /* Assemble the sequence descriptor: */
    int8_t estep = (int8_t)(-expSub);
    int32_t skip = 0;
    int32_t size = den*(ns - 1) - 2*skip + 1;
    (*seq) = msm_seq_desc_make(id, name, FALSE, size, estep, skip);
  }

void msm_test_seq_make_derived
  ( msm_seq_desc_t *seqo,
    double_vec_t *smpa,
    double mutDev,
    double delProb,
    msm_seq_id_t idd, 
    char *named,
    msm_seq_desc_t *seqd,
    double_vec_t *smpd,
    msm_rung_vec_t *gv
  )
  { msm_rung_vec_t gvr;
    /* Make a mutated copy {smpd} of the sample vector, noting the pairing {gvr}: */
    msm_double_vec_mutate(smpa, mutDev, delProb, smpd, &gvr);
    int32_t nsd = smpd->ne; /* Number of derived seq. samples. */
    /* Compute the subsampling factor {den}: */
    int8_t estepo = seqo->estep;
    int32_t skipo = seqo->skip;
    assert((estepo >= -10) && (estepo <= 0));
    int32_t den = (1 << -estepo);
    /* Assemble the sequence descriptor: */
    bool_t revd = seqo->rev;
    int8_t estepd = estepo;
    int32_t skipd = skipo;
    int32_t sized = den*(nsd - 1) - 2*skipd + 1;
    (*seqd) = msm_seq_desc_make(idd, named, revd, sized, estepd, skipd);
    /* Scale all rungs in {gvr} by {den}: */
    int32_t i;
    for (i = 0; i < gvr.ne; i++)
      { msm_rung_t *gvi = &(gvr.e[i]);
        gvi->c[0] = gvi->c[0]*den - skipo;
        gvi->c[1] = gvi->c[1]*den - skipd;
      }
    /* Interpolate those rungs: */
    (*gv) = msm_rung_vec_interpolate(&gvr);
    /* Recycle the temp storage: */
    free(gvr.e);
  }
  
void msm_fake_initial_cands
  ( int32_t ntr, 
    msm_seq_desc_t *seq0,
    msm_seq_desc_t *seq1,
    msm_cand_vec_t *cdv,
    int32_t *ncdP
  )
  { /* Sequence lengths: */
    int32_t n0 = seq0->size;
    int32_t n1 = seq1->size;
    /* Define the max candidate length {maxlen}: */
    int32_t maxlen = (n0 > n1 ? n0 : n1);
    int32_t minlen = (maxlen + 9)/10;
    /* Generate the candidates: */
    msm_cand_vec_throw
      ( ntr, seq0, seq1,
        minlen, maxlen,
        /*atomProb*/ 1.0,
        /*diagProb*/ 0.2,
        /*skipProb*/ 0.1,
        cdv, ncdP
      );
  }

double msm_test_get_sample(msm_seq_desc_t *seq, double_vec_t *smp, msm_seq_pos_t i)
  { double f = msm_seq_desc_map_index_to_orig_seq((double)i, seq);
    double v = msm_double_vec_interpolate(smp, f);
    return v;
  }

double msm_test_sample_diffsq
  ( msm_seq_desc_t *seq0,
    double_vec_t *smp0,
    msm_seq_pos_t i0,
    msm_seq_desc_t *seq1,
    double_vec_t *smp1, 
    msm_seq_pos_t i1,
    bool_t debug
  )
  { double v0 = msm_test_get_sample(seq0, smp0, i0);
    double v1 = msm_test_get_sample(seq1, smp1, i1);
    double dv = v0 - v1;
    if (debug) { fprintf(stderr, "[ v0 = %7.4f  v1 = %7.4f dv = %7.4f ]\n", v0, v1, dv); }
    return dv*dv;
  }

double msm_test_step_score
  ( msm_seq_desc_t *seq0,
    double_vec_t *smp0,
    msm_seq_desc_t *seq1,
    double_vec_t *smp1,
    msm_rung_t *g, 
    msm_rung_t *h,
    double EQL,
    double DIF,
    double BRK,
    double SKP
  )
  { 
    bool_t debug = FALSE;
    
    /* Get the indices into each sequence: */
    int32_t ig0 = g->c[0], ig1 = g->c[1];
    int32_t ih0 = h->c[0], ih1 = h->c[1];
    /* See which rungs are defined: */
    bool_t undg = msm_rung_is_none(g);
    bool_t undh = msm_rung_is_none(h);
    /* Compute a measure {Sm} of how much the step deviates from perfection: */
    double Sm;
    int32_t d0, d1;
    if (undg || undh)
      { d0 = d1 = 0; Sm = 0; }
    else 
      { d0 = ih0 - ig0; 
        d1 = ih1 - ig1;
        if ((d0 <= 0) || (d1 <= 0))
          { fprintf(stderr, " (%d %d) --> (%d %d)", ig0, ig1, ih0, ih1);
            demand (FALSE, "bad step");
          }
        Sm = 
          BRK*((d0 != 1) || (d1 != 1) ? 1.0 : 0.0) + 
          SKP*(double)(d0 + d1 - 2);
      }
    /* Compute a measure {Sd} of the difference between the paired samples: */
    double Sd;
    if (undg && undh)
      { Sd = 0; }
    else
      { /* Compute the sample diffs squared {ds0,ds1} at rungs {g,h}: */
        double Sdg = (undg ? 0 : msm_test_sample_diffsq(seq0, smp0, ig0, seq1, smp1, ig1, FALSE));
        double Sdh = (undh ? 0 : msm_test_sample_diffsq(seq0, smp0, ih0, seq1, smp1, ih1, FALSE));
        double Sdm = (Sdg + Sdh)/2;
        Sd = EQL*(1 - Sdm) + DIF*Sdm;
      }
    /* The score is the sum of the two parts: */
    double S = Sm + Sd;
    if (debug)
      { fprintf(stderr, "  step ");
        fprintf(stderr, "(%4d %4d)", (undg ? -1 : ig0), (undg ? -1 : ig1));
        if (! undg) 
          { double vg0 = msm_test_get_sample(seq0, smp0, ig0);
            double vg1 = msm_test_get_sample(seq1, smp1, ig1);
            double dgsq = msm_test_sample_diffsq(seq0, smp0, ig0, seq1, smp1, ig1, FALSE);
            fprintf(stderr, " = (%7.4f %7.4f) : %7.4f", vg0, vg1, dgsq);
          }
        else
          { fprintf(stderr, "%*s", 30, ""); }
        fprintf(stderr, " --> ");
        fprintf(stderr, "(%4d %4d)", (undh ? -1 : ih0), (undh ? -1 : ih1));
        if (! undh) 
          { double vh0 = (undh ? 0.0 : msm_test_get_sample(seq0, smp0, ih0));
            double vh1 = (undh ? 0.0 : msm_test_get_sample(seq1, smp1, ih1));
            double dhsq = msm_test_sample_diffsq(seq0, smp0, ih0, seq1, smp1, ih1, FALSE);
            fprintf(stderr, " = (%7.4f %7.4f) : %7.4f", vh0, vh1, dhsq);
          }
        else
          { fprintf(stderr, "%*s", 30, ""); }
        if (! (undg || undh))
          { fprintf(stderr, " d = (%4d %4d)", d0, d1); }
        else
          { fprintf(stderr, "%*s", 16, ""); }
        fprintf(stderr, " Sm = %12.5f Sd = %12.5f S = %12.5f", Sm, Sd, S);
        fprintf(stderr, "\n");
      }
    return S;
  }

void msm_show_tableau
  ( char *outName,
    int32_t ic, /* Candidate index. */
    int32_t ir, /* Refinement iteration index. */
    msm_cand_t *cd, /* Refined candidate. */
    msm_dyn_tableau_t *tb  /* Dynamic programming tableau. */
  )
  {
    msm_seq_desc_t *seq0 = &(cd->seq[0]);
    msm_seq_desc_t *seq1 = &(cd->seq[1]);
    char *tag = NULL;
    asprintf(&tag, "-%05d-%d", ic, ir);
    int32_t ng = msm_pairing_num_rungs(cd->pr);
    msm_rung_t gopt = msm_pairing_get_rung(cd->pr, ng-1);
    bool_t scale = FALSE;
    bool_t colorize = TRUE;
    msm_image_dyn_tableau_write_named(seq0->size, seq1->size, tb, gopt, scale, colorize, outName, tag);
    free(tag);
  }

msm_options_t *msm_get_options(int32_t argc, char**argv)
  { 
    msm_options_t *o = (msm_options_t *)notnull(malloc(sizeof(msm_options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_get_keyword(pp, "-seqLength");
    o->seqLength = (int32_t)argparser_get_next_int(pp, 1, INT32_MAX);
    
    argparser_get_keyword(pp, "-expSub");
    o->expSub = (int32_t)argparser_get_next_int(pp, 0, 10);
    
    argparser_get_keyword(pp, "-mutDev");
    o->mutDev = (int32_t)argparser_get_next_double(pp, 0.0, 1.0);
    
    argparser_get_keyword(pp, "-delProb");
    o->delProb = argparser_get_next_double(pp, 0.0, 1.0 - o->mutDev);
    
    if (argparser_keyword_present(pp, "-delta"))
      { o->delta = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->delta = 3; }
    
    if (argparser_keyword_present(pp, "-kappa"))
      { o->kappa = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->kappa = 6; }
    
    if (argparser_keyword_present(pp, "-expand"))
      { o->expand = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->expand = 0; }

    if (argparser_keyword_present(pp, "-shrink"))
      { o->shrink = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->shrink = 0; }

    if (argparser_keyword_present(pp, "-maxUnp"))
      { o->maxUnp = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->maxUnp = 6; }
    
    argparser_get_keyword(pp, "-nCands");
    o->nCands = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX);
    
    if (argparser_keyword_present(pp, "-repeat"))
      { o->repeat = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->repeat = 1; }
    
    argparser_skip_parsed(pp);
    
    o->outName = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }

#define PROG_NAME "test_dnae_filter"
#define PROG_DESC "test of DNA signal filtering routines (also spectrum plots)"
#define PROG_VERS "1.0"

/* Last edited on 2023-11-26 06:53:10 by stolfi */

#define test_dnae_filter_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -maxLevel {MAX_LEVEL} \\\n" \
  "  -initFilter " wt_table_args_parse_weights_HELP " \\\n" \
  "  -initStep {INIT_STEP} \\\n " \
  "  -incrFilter " wt_table_args_parse_weights_HELP " \\\n" \
  "  [ -plotSignals ] [ -plotSpectra ] \\\n" \
  "  [ -plotSize {PLOT_WIDTH} {PLOT_HEIGHT} ] \\\n" \
  "  [ -fontSize {FONT_SIZE} ] \\\n" \
  "  {SEQ_FILE_NAME} \\\n" \
  "  {OUT_NAME}"

#define stringify(x) #x

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program reads a" \
  " DNA/RNA sequence from the file named \"{SEQ_FILE_NAME}.bas\", converts it" \
  " to numeric format, and outputs" \
  " filtered and subsampled versions of that" \
  " sequence.  The input file should be in FASTA (\".bas\") format.\n" \
  "\n" \
  "  The program optionally generates plots of the spectra of idealized random numeric DNA" \
  " signals filtered at various scales with the given" \
  " filter coeffcients." \
  "\n" \
  "INPUT FILE FORMAT\n" \
  "  " dnae_nucleic_file_format_INFO "\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  All output files will have names starting with {OUT_NAME}.\n" \
  "\n" \
  "OPTIONS\n" \
  "\n" \
  "  -maxLevel {MAX_LEVEL}\n" \
  "    This mandatory argument specifies the" \
  " maximum level of filtering to be" \
  " applied.  Level 0 is the unfiltered" \
  " sequence, and each subsequent level is sampled at half" \
  " the frequency of the previous one.\n" \
  "\n" \
  "  -initFilter " wt_table_args_parse_weights_HELP " \n" \
  "  -incrFilter " wt_table_args_parse_weights_HELP " \n" \
  "    These mandatory arguments specify the" \
  " weights of the filter to be used at the" \
  " first filtering step (from level 0 to level 1) and at subsequent" \
  " filtering steps.  " wt_table_args_parse_weights_norm_sum_INFO "\n" \
  "\n" \
  "  -initStep {INIT_STEP} \n" \
  "    This mandatory argument specifies the initial resampling step, which must" \
  " be a power of 2.  If positive, the sequence will be downsampled after" \
  " filtering, retaining only one every {INIT_STEP} datums.  If negative, the" \
  " sequence will be interpolated with new datums spaced {INIT_STEP} times original step.\n" \
  "\n" \
  "  -plotsignals \n" \
  "    This optional argument requests plots of the filtered numeric encoding, as" \
  " Encapsulated PostScript (\".eps\") files, with each channel in a distinct" \
  " color (red, green, and blue).\n" \
  "\n" \
  "  -plotSpectra \n" \
  "    This optional argument requests plots of the spectra of idealized random" \
  " sequences as" \
  " Encapsulated PostScript (\".eps\") files.  For this computation, the sequence" \
  " length is rounded to the next power of two, and its input spectrum is assumed" \
  " to have uniform power at all frequencies.  The filtering kernel is assumed" \
  " to wrap around so that there is no loss of datums at the ends. \n" \
  "\n" \
  "  -plotSize {PLOT_WIDTH} {PLOT_HEIGHT} \n" \
  "    This optional argument specifies the width and height of the" \
  " plot (excluding scales) in millimeters.  If" \
  " omitted, the program assumes " stringify(dnae_DEFAULT_PLOT_WIDTH) "" \
  " by " stringify(dnae_DEFAULT_PLOT_HEIGHT) " mm.\n" \
  "\n" \
  "  -fontSize {FONT_SIZE} \n" \
  "    This optional argument specifies the font size in points (pt).  If" \
  " omitted, the program assumes " stringify(dnae_DEFAULT_FONT_SZ) " pt.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  The light at the end of the tunnel.\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 2006-12-16 by J. Stolfi as {dm_test_010_filter.c}.\n" \
  "  Parts of it were created 2006-12-27 by J. Stolfi as {dm_test_012_fourier.c}.\n" \
  "MODIFICATION HISTORY\n" \
  "  2014-06-10 J. Stolfi: renamed {test_dnae_filter}. \n" \
  "  2014-07-29 J. Stolfi: added the \"-fontSize\" and \"-plotSize\" options.\n" \
  "  2014-08-05 J. Stolfi: merged {dm_test_012_fourier.c} into this code.\n" \
  "  2014-08-10 J. Stolfi: rounded sequences to power of 2 for the power spectra plots.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_dnae_filter_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <fftw3.h>

#include <argparser.h>
#include <affirm.h>
#include <vec.h>
#include <epswr.h>
#include <wt_table.h>
#include <wt_table_args_parse.h>
#include <rn.h>

#include <msm_multi.h>
#include <msm_ps_tools.h>
#include <msm_test_tools.h>

#include <dnae_nucleic.h>
#include <dnae_seq.h>
#include <dnae_spectrum.h>
#include <dnae_test_tools.h>
#include <dnae_seq_multi.h>

#define dnae_DEFAULT_PLOT_WIDTH 180.0
#define dnae_DEFAULT_PLOT_HEIGHT 60.0
#define dnae_DEFAULT_FONT_SIZE 12.0

typedef struct options_t 
  { char *seqFileName;   /* Name of file with DNA/RNA sequence (minus ".bas" ext). */
    int32_t maxLevel;        /* Maximum filtering level. */
    double_vec_t w0;     /* Filter weights for first stage (level 0 to level 1). */
    int8_t initStep_exp; /* Exponent of 2 of resampling step from level 0 to level 1. */ 
    double_vec_t w1;     /* Filter weights for subsequent stages. */
    bool_t plotSignals;  /* TRUE to plot signals. */
    bool_t plotSpectra;  /* True to plot power spectra. */
    double fontSize;     /* The font size to use in the plots (pt). */
    double plotSize_h;   /* Width of plot area (mm). */
    double plotSize_v;   /* Height of plot area (mm). */
    char *outName;       /* Output file name prefix (minus extensions). */
  } options_t;
  
int32_t main(int32_t argc, char**argv);

options_t *dnae_get_options(int32_t argc, char**argv);
  /* Parses the command line options, packs 
    them into a {options_t} record. */

void dnae_generate_and_plot_spectra(options_t *o, int32_t N);
  /* Generates and plots ideal spectra of random DNA signals
    and their filtered/sampled versions, assuming {N}
    datums in the basic sequence. */

void dnae_compute_random_dna_spectrum(int32_t N, double randP[]);
  /* Stores in {randP[0..N/2]} the theoretical power spectrum of a
    random DNA sequence. Assumes that the bases are uniformly
    distributed and independent. Also assumes that they are converted
    to numeric 3-vectors {(±1,±1,±1)} as described under
    {dnae_datum_decoded_from_nucleic_char}. The number {N} must be even.
    
    The spectrum is the sum of the spectra of all three channels */

void dnae_compute_filter_spectrum(int32_t nw, double w[ ], int32_t N, double wP[]);
  /* Stores in {wP[0..N/2]} the power spectrum of the sequence {w[0..nw-1]},
    zero-extended to {N} samples.  Requires {nw <= N}. */

void dnae_dump_dft(int32_t N, char *xN, fftw_complex *xP, char *yN, fftw_complex *yP);
  /* Writes the complex vectors {xP[0..N-1]} and {yP[0..N-1]}
    (typically, a complex discrete signal and its DFT), side by side.
    The parameters {xN,yN} should be the names of the three vectors.
    If either vector is NULL, it is omitted from the printout. */

void dnae_dump_spectra(int32_t fMax, char *xN, double *xP, char *yN, double *yP, char *zN, double *zP);
  /* Writes the power spectra {xP[0..fMax]}, {yP[0..fMax]},
    {zP[0..fMax]} side by side. The parameters {xN,yN,zN} should be
    the names of the three vectors. If either vector is NULL, it is
    omitted from the printout. */

void dnae_plot_spectra
  ( int32_t fMax,
    double *xP,
    double *yP,
    double *zP,
    int32_t fSep,
    double vMax,
    bool_t skipZero,
    double hPlotSize,
    double vPlotSize,
    double fontSize,
    char *outName,
    int32_t level,
    char *tag
  );
  /* Plots up to three power spectra on the same graph, using the same 
    scales.  
    
    The X axis runs from 0 to {fMax}, the Y axis from 0 to {vMax} or
    the maximum power in the data, whichever is greater.
    
    The power spectra are assumed to be {xP[0..fMax]},
    {yP[0..fMax]}, {zP[0..fMax]}.  The zero-frequency point is 
    omitted if {skipzero} is true.  If any of {xP,yP,zP} is NULL, the 
    corresponding graph is omitted.  A vertical line is drawn at 
    frequency {fSep}.  The plot is written to file
    "{outName}-{level}{tag}.eps". */

void dnae_plot_graph(msm_ps_tools_t *dp, int32_t fMax, double xP[], bool_t skipZero, double R, double G, double B, bool_t show_dots);
  /* Plots to {dp} the graph with points {(f,xP[f])} for {f} between {0}
    and {fMax} inclusive. Note that {xP} must have {fMax+1} elements.
    The first data point is omitted if {skipZero} is true. The
    points are connecetd by a polyline with the given {R,G,B} color; if
    {show_dots} is true, also a dot is drawn at each point. */

options_t *dnae_get_options(int32_t argc, char**argv);
  /* Parses the command line options, packs 
    them into a {options_t} record. */

#define dnae_debug_fft FALSE

int32_t main(int32_t argc, char**argv)
  { 
    options_t *o = dnae_get_options(argc, argv);

    fprintf(stderr, "reading input sequence ...\n");
    dnae_seq_t s = dnae_seq_read_from_nucleic_file_named(o->seqFileName, "", ".bas", 0, "s", FALSE);
    fprintf(stderr, "sequence has %d sample datums\n", dnae_seq_num_datums(&s));
    
    fprintf(stderr, "getting the weight tables ...\n");
    char *wname0 = wt_table_make_descr(o->w0.ne, o->w0.e, "%6.4f");
    wt_table_print(stderr, "ini", o->w0.ne, o->w0.e, 0);
    char *wname1 = wt_table_make_descr(o->w1.ne, o->w1.e, "%6.4f");
    wt_table_print(stderr, "inc", o->w1.ne, o->w1.e, 0);

    fprintf(stderr, "trying to filter up to level %d ...\n", o->maxLevel);
    dnae_seq_t r[o->maxLevel + 1];
    int8_t ek0 = 0; /* No resampling going from level 0 to level 1. */
    dnae_seq_multi_filter(&s, o->maxLevel, &(o->w0), wname0, &(o->w1), wname1, ek0, r);
    
    /* Figure out max useful level {maxLevel}: */
    int32_t maxLevel = o->maxLevel;
    while ((maxLevel >= 0) && (dnae_seq_num_datums(&(r[maxLevel])) == 0)) { maxLevel--; }
    fprintf(stderr, "maximum non-empty level is %d\n", maxLevel);
    
    fprintf(stderr, "writing and plotting the filtered sequences ...\n");
    dnae_test_tools_seq_multi_write_and_plot_named
      ( r, maxLevel, "test seq", o->outName, "", 
        o->plotSignals, o->plotSize_h, o->plotSize_v, o->fontSize, -1
      );

    if (o->plotSpectra) 
      { int32_t N = 1;
        while (N < s.sd.size) { N <<= 1; }
        dnae_generate_and_plot_spectra(o, N);
      }
  
    /* ...and we are done: */
    fprintf(stderr, "done.\n");
    return 0;
  }
      
void dnae_generate_and_plot_spectra(options_t *o, int32_t N)
  {
    /* Convenience: */
    double hSz = o->plotSize_h;
    double vSz = o->plotSize_v;
    double fSz = o->fontSize;
    int32_t ek0 = o->initStep_exp;
    
    /* Number of datums in original signal: */
    int32_t origN = N;
    
    /* Maximum frequency in original signal: */
    int32_t origF = origN/2; 
    
    /* Nominal maximum power in random DNA spectrum: */
    double origW = 6.0*dnae_NUCLEIC_RAW_SCALE*dnae_NUCLEIC_RAW_SCALE;

    /* Power spectrum of original sequence (simulated): */
    fprintf(stderr, "computing the mean power spectrum of a random periodic DNA sequence with %d datums...\n", N);
    double *origP = notnull(malloc((origF+1)*sizeof(double)), "no mem"); /* Non-aliased part of {fP}. */
    dnae_compute_random_dna_spectrum(origN, origP);
    dnae_plot_spectra(origF, origP, NULL, NULL, 0, origW, FALSE, hSz, vSz, fSz, o->outName, 0, "-db");

    /* Maximum samples any filtered sequence (accounting for subsampling at level 1): */
    int32_t maxN = origN*(ek0 < 0 ? (1 << (-ek0)) : 1); 
    
    /* Maximum frequency in any power spectra: */
    int32_t maxF = maxN/2;
    
    /* Nominal maximum power in filter spectra: */
    double maxW = 1.0;
    
    fprintf(stderr, "computing the power spectra of the two filters at max resolution ...\n");
    /* wt_table_print(stderr, "initial", o->w0.ne, o->w0.e, 0); */
    /* wt_table_print(stderr, "incremental", o->w1.ne, o->w1.e, 0); */
    double *w0P = notnull(malloc((maxF+1)*sizeof(double)), "no mem"); /* Initial filter. */;
    double *w1P = notnull(malloc((maxF+1)*sizeof(double)), "no mem"); /* Incremental filter. */;
    dnae_compute_filter_spectrum(o->w0.ne, o->w0.e, maxN, w0P);
    dnae_compute_filter_spectrum(o->w1.ne, o->w1.e, maxN, w1P);
    if (dnae_debug_fft) { dnae_dump_spectra(maxF, "w0P ", w0P, "w1P ", w1P, NULL, NULL); }
    dnae_plot_spectra(maxF, w0P, NULL, NULL, 0, maxW, FALSE, hSz, vSz, fSz, o->outName, 0, "-f0");
    dnae_plot_spectra(maxF, w1P, NULL, NULL, 0, maxW, FALSE, hSz, vSz, fSz, o->outName, 0, "-f1");
    

    fprintf(stderr, "running Fourier analysis up to level %d ...\n", o->maxLevel);
    
    /* !!! Should account for quantization error too. !!! */
    
    /* !!! Should account for sequence shortening. !!! */
    
    int32_t f;
    
    /* Initialize spectra for original sequence (level 0): */
    int32_t curL = 0;     /* Current filtering level. */
    int32_t curN = origN;  /* Number of samples in current signal. */
    int32_t curF = origF;  /* Max freq in {vP,aP}. */
    
    /* Spectrum of filtered DNA signal at max resolution, assuming no resampling: */
    double *fP = rn_alloc((maxF+1));
    for (f = 0; f <= maxF; f++) { fP[f] = (f <= origF ? origP[f] : 0.0); }

    /* Spectra of components of the current signal (filtered and resampled): */
    double *vP = notnull(malloc((curF+1)*sizeof(double)), "no mem"); /* Non-aliased part. */
    double *aP = notnull(malloc((curF+1)*sizeof(double)), "no mem"); /* Aliased part. */
    for (f = 0; f <= curF; f++) { vP[f] = fP[f]; aP[f] = 0; }
    
    /* Effective filter kernel, at max sample count: */
    double *wft = rn_alloc(maxN);
    
    /* Spectrum of effective filter kernel, at max resolution: */
    double *wfP = rn_alloc(maxF+1);
    
    /* Simulate filtering and plot spectra: */
    while (TRUE)
      { assert(curF <= maxF);
        
        /* Here {fP[0..maxF]} is the power spectrum of the original signal with all the filtering
          steps applied, but without any downsampling.  Assumes that any interpolation 
          is perfect so that it preserves the original power spectrum at low frequencies and
          introduces no high-frequency terms. Also {vP[0..curF]} and {aP[0..curF]}
          are power spectra of two parts of the current signal (filtered and resampled):
          {vP} is all the Fourier terms that have not been aliased and {aP}
          are the aliased terms, mapped to the lower equivalent frequency. */
        
        /* Plot spectrum of filtered, non-sampled signals: */
        /* Draw vertical line at the cutoff freq {curF} of this level. */
        fprintf(stderr, "plotting full spectrum of filtered signal\n");
        if (dnae_debug_fft) { dnae_dump_spectra(maxF, "fP", fP, NULL, NULL, NULL, NULL); }
        dnae_plot_spectra(maxF, fP, NULL, NULL, curF, maxW, TRUE, hSz, vSz, fSz, o->outName, curL, "-fo");

        /* Compute the rms frequency and compare to the current cutoff frequency: */
        double sum_Pf2 = 0;
        double sum_P = 0;
        for (f = 0; f <= maxF; f++) { sum_P += fP[f]; sum_Pf2 += fP[f]*f*f; }
        double rmsF = sqrt(sum_Pf2/sum_P);
        fprintf(stderr, "rms frequency at full resolution = %8.4f (%8.6f*curF)\n", rmsF, rmsF/curF);
        
        /* Plot spectra of non-aliased and aliased filtered and sampled signals: */
        fprintf(stderr, "plotting spectra of non-aliased and aliased parts of filtered sampled signals:\n");
        if (dnae_debug_fft) { dnae_dump_spectra(curF, "vP", vP, "aP", aP,   NULL, NULL); }
        dnae_plot_spectra(curF, vP, aP,   NULL, 0,    maxW, TRUE, hSz, vSz, fSz, o->outName, curL, "-fs");

        /* Compute the current signal and noise powers and N/S ratio: */
        double vPt = 0; /* Total power in non-aliased freqencies. */
        double aPt = 0; /* Total power in aliased frequencies. */
        int32_t g;
        for (g = 0; g <= curF; g++) { vPt += vP[g]; aPt += aP[g]; }
        fprintf
          ( stderr,
            "total power S+N = %8.5f  valid S = %8.5f  aliased N = %8.5f  N/(S+N) = %8.6f\n",
            vPt+aPt, vPt, aPt, aPt/(vPt+aPt)
          );

        /* Can we do another stage of filtering? */
        if (curL >= o->maxLevel) { break; }
        if (o->w1.ne > curN) { fprintf(stderr, "!! seq too short\n"); break; }

        /* Compute number of samples {nextN} and max freq {nextF} after next resampling: */
        int32_t nextN;
        if (curL == 0)
          { if (ek0 < 0) 
              { /* Interpolation at level 0 --> level 1: */
                nextN = curN*(1 << (-ek0));
              }
            else if (ek0 > 0)
              { /* Downsampling at level 0 --> level 1: */
                int32_t step = (1 << ek0);
                if (curN % step != 0) { fprintf(stderr, "!! cannot evenly downsample\n"); break; }
                nextN = curN/step;
              }
            else
              { nextN = curN; }
          }
        else
          { assert(curN % 2 == 0);
            nextN = curN/2;
          }
        int32_t nextF = nextN/2; /* Max freq after next resampling. */

        /* APPLY FILTERING: */
        
        /* Get the proper filter {wt} to use: */
        double_vec_t *wt = (curL == 0 ? &(o->w0) : &(o->w1)); /* Filter weight kernel to use. */
        
        /* Build effective kernel {wft[0..maxN-1]} of filter in max resolution: */
        assert(maxN % curN == 0);
        int32_t maxS = maxN/curN; /* Step of current signal relative to max resolution signal: */
        assert(wt->ne*maxS <= maxN);
        int32_t t;
        for (t = 0; t < maxN; t++) { wft[t] = 0; }
        for (t = 0; t < wt->ne; t++) { wft[t*maxS] = wt->e[t]; }
        /* Compute power spectrum of incremental filter in max resolution: */
        fprintf(stderr, "computing full-res spectrum of incremental filter for level %d --> %d:\n", curL, curL+1);
        dnae_compute_filter_spectrum(maxN, wft, maxN, wfP);
        if (dnae_debug_fft) { dnae_dump_spectra(maxF, "wfP", wfP, NULL, NULL, NULL, NULL); }

        /* Apply incremental filter to unsampled signal {fP}: */
        for (f = 0; f <= maxF; f++) { fP[f] *= wfP[f]; }
    
        /* Compute power spectrum of incremental filter in current resolution: */
        double *wcP = rn_alloc(curF+1);
        fprintf(stderr, "computing current-res spectrum of incremental filter for level %d --> %d:\n", curL, curL+1);
        dnae_compute_filter_spectrum(wt->ne, wt->e, curN, wcP);
        dnae_plot_spectra(curF, wcP, NULL, NULL, nextF, maxW, TRUE, hSz, vSz, fSz, o->outName, curL+1, "-fr");
        if (dnae_debug_fft) { dnae_dump_spectra(curF, "wcP", wcP, NULL, NULL, NULL, NULL); }

        /* Apply incremental filter to {vP,aP}: */
        for (g = 0; g <= curF; g++) { vP[g] *= wcP[g]; aP[g] *= wcP[g]; }
        free(wcP); 

        /* APPLY RESAMPLING: */
        
        /* Allocate new storage for valid and aliased spectra after next resampling: */
        double *nextvP = notnull(malloc((nextF+1)*sizeof(double)), "no mem");
        double *nextaP = notnull(malloc((nextF+1)*sizeof(double)), "no mem");
        
        /* Simulate resampling of {vP,aP} to {nextvP,nextaP}: */
        int32_t h;
        for (h = 0; h <= nextF; h++) { nextvP[h] = nextaP[h] = 0; }
        if (nextN >= curN)
          { /* No downsampling, hence no aliasing: */
            for (g = 0; g <= curF; g++) { nextvP[g] = vP[g]; nextaP[g] = aP[g]; }
          }
        else
          { assert(curN % nextN == 0);
            int32_t nextS = curN/nextN; /* Step of next signal relative to current one. */
            for (g = 0; g <= curF; g++) 
              { /* Compute the frequency {h} that {g} becomes when 2-sampled: */
                h = nextS*g;
                if (h > curF) h = h - nextS*curF;
                if (h < 0) { h = -h; }
                assert(h % nextS == 0);
                h /= nextS;
                assert((h >= 0) && (h <= nextF));
                /* Add power {vP[g]}  to {nextvP[g]} or {nextaP[g]} as appropriate: */
                if (g == h) 
                  { /* Freq {g} gets preserved: */      nextvP[h] += vP[g]; }
                else
                  { /* Freq {g} gets aliased to {h}: */ nextaP[h] += vP[g]; }
                /* Add the power that was aliased in previous levels: */
                nextaP[h] += aP[g];
              }
          }
          
        /* Prepare for next iteration: */
        free(vP); vP = nextvP;
        free(aP); aP = nextaP;
        curL++;
        curN = nextN;
        curF = nextF;
      }

    /* Cleanup: */
    free(origP);
    free(w0P);
    free(w1P);
    free(fP);
    free(vP);
    free(aP);
    free(wft);
    free(wfP);
  }

void dnae_dump_spectra(int32_t fMax, char *xN, double *xP, char *yN, double *yP, char *zN, double *zP)
  { int32_t f;
    for (f = 0; f <= fMax; f++)
      { if (xP != NULL) { fprintf(stderr, " %s[%4d] = %8.6f", xN, f, xP[f]); }
        if (yP != NULL) { fprintf(stderr, " %s[%4d] = %8.6f", yN, f, yP[f]); }
        if (zP != NULL) { fprintf(stderr, " %s[%4d] = %8.6f", zN, f, zP[f]); }
        fprintf(stderr, "\n");
      }
  }

void dnae_dump_dft(int32_t N, char *xN, fftw_complex *xP, char *yN, fftw_complex *yP)
  { int32_t i;
    for (i = 0; i < N; i++)
      { if (xP != NULL)
          { fprintf(stderr, "  %s[%4d] = ( %+9.6f %+9.6f )", xN, i, xP[i][0], xP[i][1]); }
        if (yP != NULL) 
          { fprintf(stderr, "  %s[%4d] = ( %+9.6f %+9.6f )", yN, i, yP[i][0], yP[i][1]); }
        fprintf(stderr, "\n"); 
      }
  }

void dnae_compute_random_dna_spectrum(int32_t N, double randP[])
  { /* Get and check some parameters: */
    demand(N % 2 == 0, "N must be even");
    int32_t fMax = N/2; /* Maximum frequency n power spectrum. */
    int32_t nc = dnae_CHANNELS; /* Number of channels in datums. */
    assert(nc == 3); /* The analysis below is specific for {nc=3}. */
    double uP = 1.0; /* Expected power per channel and per freq {0..N-1} */
    /* In each channel, the samples are {-1} or {+1} with equal prob. */
    /* The expected value of the mean squared is {1/N}: */
    randP[0] = nc*uP;
    /* For frequency {fMax}, the expected power is {1/N} too: */
    randP[fMax] = nc*uP;
    /* For other freqs {1..fMax-1}, the expected power is {2/N}  too: */
    int32_t f;
    for (f = 1; f < fMax; f++) { randP[f] = 6*uP; }
  }

void dnae_compute_filter_spectrum(int32_t nw, double w[ ], int32_t N, double wP[])
  { /* Get and check some parameters: */
    demand(N % 2 == 0, "N must be even");
    int32_t fMax = N/2; /* Maximum frequency in power spectrum. */
    demand(nw <= N, "nw must not exceed N");
    
    /* Allocate vectors {in,ot} for {fftw3.h} routines: */
    fftw_complex *in = fftw_malloc(sizeof(fftw_complex)*N); /* Input vector. */
    fftw_complex *ot = fftw_malloc(sizeof(fftw_complex)*N); /* Output vector. */

    /* Precompute parameters for {fftw3.h} routines: */
    fftw_plan plan = fftw_plan_dft_1d(N, in, ot, FFTW_FORWARD, FFTW_ESTIMATE);

    /* Store weights {w} in input vector, as complex numbers: */
    int32_t i;
    for (i = 0; i < N; i++)
      { in[i][0] = (i < nw ? w[i] : 0); 
        in[i][1] = 0; 
      }

    /* Compute Fourier transform: */
    fftw_execute(plan);
      
    if (dnae_debug_fft)
      { fprintf(stderr, "input and output of Fouier transform:\n");
        dnae_dump_dft(N, "in", in, "ot", ot);
      }

    /* Compute the power spectrum: */
    double norm = 1/((double)N); /* Power normalization factor. */
    int32_t f;
    for (f = 0; f <= fMax; f++)
      { double P = ot[f][0]*ot[f][0] + ot[f][1]*ot[f][1];
        if ((f > 0) && (f < fMax))
          { int32_t g = N-f;
            P += ot[g][0]*ot[g][0] + ot[g][1]*ot[g][1];
          }
        wP[f] = norm*P;
      }
      
    if (dnae_debug_fft)
      { fprintf(stderr, "raw power spectrum:\n");
        for (f = 0; f <= fMax; f++)
          { fprintf(stderr, " wP[%4d] = %8.6f\n", f, wP[f]); }
      }

    /* Normalize power spectrum so that the max elem is 1: */
    double vmax = -INF;
    for (f = 0; f <= fMax; f++) { if (wP[f] > vmax) { vmax = wP[f]; } }
    for (f = 0; f <= fMax; f++) { wP[f] /= vmax; }
      
    /* Cleanup: */
    fftw_destroy_plan(plan);
    fftw_free(in); 
    fftw_free(ot);
  }

void dnae_plot_spectra
  ( int32_t fMax,
    double *xP,
    double *yP,
    double *zP,
    int32_t fSep,
    double vMax,
    bool_t skipZero,
    double hPlotSize,
    double vPlotSize,
    double fontSize,
    char *outName,
    int32_t level,
    char *tag
  )
  { /* Colors for curves (all with brightness = 0.30): */
    double R[3] = { 1.000, 0.000, 0.000 };
    double G[3] = { 0.000, 0.500, 0.330 };
    double B[3] = { 0.000, 0.000, 1.000 };
    
    /* Color for vertical bar: */
    double C[3] = { 0.800,0.500,0.000 };
    
    /* Max label chars in each scale: */
    int32_t maxXLabChars = (int32_t)ceil(log(fmax(fMax,9))/M_LN10);
    int32_t maxYLabChars = 5; /* E.g. "+1.00". */
    
    /* Figure dimensions in mm: */
    double mm_per_pt = 25.4/72.0;  /* Size of a point in mm. */
    double mrgL =  0.75*maxYLabChars*fontSize*mm_per_pt;  /* Width of Y scale and ticks. */
    double mrgR =  0.40*maxXLabChars*fontSize*mm_per_pt;  /* Possible right overflow of X scale. */
    double mrgB =  1.25*fontSize*mm_per_pt;  /* Height of X scale and ticks. */
    double mrgT =  0.25*fontSize*mm_per_pt;  /* Possible top overflow of Y scale. */
    double hTotSize = hPlotSize + mrgL + mrgR; /* Tot width of plottable area (incl scales). */
    double vTotSize = vPlotSize + mrgB + mrgT; /* Tot height of plottable area (incl scales). */
    double mrgE = 0.5; /* Extra margin around plottable area. */

    /* Create plotter object {dp} and get its Postscript stream {eps}: */
    char *levTag = NULL;
    char *levTag = jsprintf("-%02d%s", level, tag);
    msm_ps_tools_t *dp = 
      msm_ps_tools_new
        ( NULL, outName, levTag, hTotSize, vTotSize, 
          fontSize, maxXLabChars, maxYLabChars, mrgE
        );
    epswr_figure_t *eps = msm_ps_tools_get_eps_figure(dp);
    
    /* Choose nominal X range {[xMin_xMax]}: */
    double xSkosh = 0.0; /* Was  0.05*fMax;; */
    double xMin = 0.0; /* Was 0.0 - xSkosh; */
    double xMax = fMax + xSkosh;
    
    /* Enlarge {vmax} if needed to enclose all data points: */
    int32_t f;
    for (f = 0; f <= fMax; f++)
      { if ((xP != NULL) && (xP[f] > vMax)) { vMax = xP[f]; } 
        if ((yP != NULL) && (yP[f] > vMax)) { vMax = yP[f]; } 
        if ((zP != NULL) && (zP[f] > vMax)) { vMax = zP[f]; } 
      }

    /* Compute nominal Y range {[yMin_yMax]}: */
    double ySkosh = 0.10*vMax;
    double yMax = vMax + ySkosh;
    double yMin = 0.0; /* Was 0.0 - ySkosh; */
    
    /* Set Graph reference window, with some skosh: */
    msm_ps_tools_set_graph_ref_window(dp, xMin, xMax, yMin, yMax);
    
    /* Shrink the Epswr window to leave space for tick labels: */
    msm_ps_tools_shrink_epswr_ref_window(dp, mrgL, mrgR, mrgB, mrgT);
    
    /* Decide whether to plot dots at individual samples: */
    double hStep = msm_ps_tools_map_x(dp, 1) - msm_ps_tools_map_x(dp, 0);
    bool_t show_dots = (hStep >= 1.5);  /* If spaced at least 1.5 mm */

    /* Draw axes, ticks, labels, etc: */
    char *font = "Times-Roman";
    double tickSize = 0.30*fontSize*mm_per_pt;  /* In mm. */
    epswr_set_label_font(eps, font, fontSize);
    /* msm_ps_tools_draw_ref_axis(dp, epswr_axis_HOR, 0.5,0.5,0.5); */
    /* msm_ps_tools_draw_ref_axis(dp, epswr_axis_VER, 0.5,0.5,0.5); */
    msm_ps_tools_draw_ref_frame(dp, 0.5,0.5,0.5);
    msm_ps_tools_draw_scale(dp, epswr_axis_HOR, yMin, tickSize,0.0, 5.0,1.0, 0.5,0.5,0.5, "%.0f",  1.2, 12.0,0.0);
    msm_ps_tools_draw_scale(dp, epswr_axis_VER, xMin, tickSize,0.0, 3.0,0.0, 0.5,0.5,0.5, "%+.2f", 1.2,  4.0,0.0);
    
    if ((fSep > 0) && (fSep < fMax))
      { /* Draw a vertical line at {fSep}: */
        epswr_set_pen(eps, C[0],C[1],C[2], 0.25, 0.0, 0.0);
        msm_ps_tools_draw_segment(dp, (double)fSep, yMin, (double)fSep, yMax);
      }
    
    /* Draw polylines: */
    if (xP != NULL) { dnae_plot_graph(dp, fMax, xP, skipZero, R[0], G[0], B[0], show_dots); }
    if (yP != NULL) { dnae_plot_graph(dp, fMax, yP, skipZero, R[1], G[1], B[1], show_dots); }
    if (zP != NULL) { dnae_plot_graph(dp, fMax, zP, skipZero, R[2], G[2], B[2], show_dots); }

    msm_ps_tools_close(dp);
    free(levTag);
  }
  
void dnae_plot_graph(msm_ps_tools_t *dp, int32_t fMax, double xP[], bool_t skipZero, double R, double G, double B, bool_t show_dots)
  { assert(fMax >= 0);
    int32_t skip = (skipZero ? 1 : 0); /* How many data points to skip from start. */
    epswr_figure_t *eps = msm_ps_tools_get_eps_figure(dp);
    epswr_set_pen(eps, R,G,B, 0.25, 0.0, 0.0);
    msm_ps_tools_draw_y_polyline(dp, skip, fMax, xP+skip, fMax+1-skip);
    if (show_dots) 
      { epswr_set_fill_color(eps, R,G,B);
        msm_ps_tools_draw_y_dots(dp, skip, fMax, xP+skip, fMax+1-skip, 0.5, TRUE, FALSE);
      }
  }

options_t* dnae_get_options(int32_t argc, char**argv)
  { options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    argparser_get_keyword(pp, "-initFilter");
    o->w0 = wt_table_args_parse_weights(pp, TRUE);

    argparser_get_keyword(pp, "-incrFilter");
    o->w1 = wt_table_args_parse_weights(pp, TRUE);

    argparser_get_keyword(pp, "-maxLevel");
    o->maxLevel = (int32_t)argparser_get_next_int(pp, 0, 20);

    if (argparser_keyword_present(pp, "-initStep"))
      { double st = argparser_get_next_double(pp, 1.0/128.0, 128.0);
        int32_t ek = 0;
        while (st < 1) { ek--; st = st*2; }
        while (st > 1) { ek++; st = st/2; }
        if (st != 1.0) { argparser_error(pp, "initial sampling step must be a power of 2."); }
        assert(abs(ek) <= msm_seq_desc_estep_MAX);
        o->initStep_exp = (int8_t)ek;
      }
    else
      { /* Subsample with step 1: */
        o->initStep_exp = 0;
      }

    o->plotSignals = argparser_keyword_present(pp, "-plotSignals");

    o->plotSpectra = argparser_keyword_present(pp, "-plotSpectra");

    if (argparser_keyword_present(pp, "-plotSize"))
      { o->plotSize_h = argparser_get_next_double(pp, 20.0, 1000.0);
        o->plotSize_v = argparser_get_next_double(pp, 5.0, 1000.0);
      }
    else
      { o->plotSize_h = dnae_DEFAULT_PLOT_WIDTH;
        o->plotSize_v = dnae_DEFAULT_PLOT_HEIGHT;
      }

    if (argparser_keyword_present(pp, "-fontSize"))
      { o->fontSize = argparser_get_next_double(pp, 4.0, 240.0); }
    else
      { o->fontSize = dnae_DEFAULT_FONT_SIZE; }

    argparser_skip_parsed(pp);
    
    o->seqFileName = argparser_get_next(pp);
    
    o->outName = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }

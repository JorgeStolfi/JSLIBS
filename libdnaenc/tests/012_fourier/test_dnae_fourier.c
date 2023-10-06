#define PROG_NAME "test_dnae_fourier"
#define PROG_DESC "Fourier analysis of DNA signal filtering routines"
#define PROG_VERS "1.0"

/* Last edited on 2023-10-01 19:48:30 by stolfi */

#define test_dnae_fourier_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -maxLevel {MAX_LEVEL} \\\n" \
  "  -numPoints {N_POINTS} \\\n" \
  "  -initFilter " wt_table_args_HELP " \\\n" \
  "  -incrFilter " wt_table_args_HELP " \\\n" \
  "  [ -plotSize {PLOT_WIDTH} {PLOT_HEIGHT} ] \\\n" \
  "  [ -fontSize {FONT_SIZE} ] \\\n" \
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
  "  This program analyzes the spectrum of numeric DNA" \
  " signals filtered at various scales with given" \
  " filter coeffcients." \
  " .\n" \
  "\n" \
  "INPUT FILE FORMAT\n" \
  dnae_nucleic_file_format_INFO "\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  All output files will have names starting with {OUT_NAME}.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -maxLevel {MAX_LEVEL}\n" \
  "    This mandatory argument specifies the" \
  " maximum level of filtering to be" \
  " applied.  Level 0 is the unfiltered" \
  " sequence, and each subsequent level is sampled at half" \
  " the frequency of the previous one.\n" \
  "\n" \
  "  -numPoints {N_POINTS}\n" \
  "    This mandatory argument specifies the" \
  " number of data points (nucleotide bases)" \
  " to use at the lowest (unfiltered) scale.\n" \
  "\n" \
  "  -initFilter " wt_table_args_HELP " \n" \
  "  -incrFilter " wt_table_args_HELP " \n" \
  "    These mandatory arguments specify the" \
  " weights of the filter to be used at the" \
  " first filtering step and at subsequent" \
  " filtering steps.  " wt_table_args_norm_sum_INFO "\n" \
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
  "  dm_match(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 27/dec/2006 by J. Stolfi as {dm_test_012_fourier}.\n" \
  "MODIFICATION HISTORY\n" \
  "  2014-06-10 J. Stolfi: renamed {test_dnae_fourier}.\n" \
  "  2014-07-29 J. Stolfi: added the \"-fontSize\" and \"-plotSize\" options.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_dnae_fourier_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <fftw3.h>

#include <argparser.h>
#include <affirm.h>
#include <vec.h>
#include <wt_table.h>
#include <rn.h>

#include <msm_multi.h>
#include <msm_ps_tools.h>

#include <dnae_test_tools.h>
#include <dnae_nucleic.h>
#include <dnae_seq.h>

#define dnae_DEFAULT_PLOT_WIDTH 180.0
#define dnae_DEFAULT_PLOT_HEIGHT 60.0
#define dnae_DEFAULT_FONT_SIZE 12.0

typedef struct options_t 
  { int32_t maxLevel;       /* Maximum filtering level. */
    int32_t numPoints;      /* Number of bases in test sequence. */
    double_vec_t w0;    /* Filter weights for first stage (level 0 to level 1). */
    double_vec_t w1;    /* Filter weights for subsequent stages. */
    double fontSize;    /* The font size to use in the plots (pt). */
    double plotSize_h;  /* Width of plot area (mm). */
    double plotSize_v;  /* Height of plot area (mm). */
    char *outName;      /* Output file name prefix (minus extensions). */
  } options_t;
  
int32_t main(int32_t argc, char**argv);

void dnae_compute_random_dna_spectrum(int32_t N, double dnaP[]);
  /* Stores in {dnaP[0..N/2]} the theoretical power spectrum of a
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
  ( double *xP,
    double *yP,
    double *zP,
    int32_t fMax,
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

    /* Number of points at base level: */
    int32_t N = o->numPoints;
    
    /* Convenience: */
    double hSz = o->plotSize_h;
    double vSz = o->plotSize_v;
    double fSz = o->fontSize;
    
    fprintf(stderr, "reading input sequence ...\n");
    /* dnae_seq_t s = dnae_seq_read_from_nucleic_file(0, "s", o->seqName); */

    /* Maximum frequency in power spectra: */
    int32_t fMax = N/2;
    
    /* Nominal maximum power in random DNA spectrum: */
    double dnaMax = 6.0*dnae_NUCLEIC_RAW_SCALE*dnae_NUCLEIC_RAW_SCALE;
    
    /* Nominal maximum power in filter spectra: */
    double wMax = 1.0;
    
    fprintf(stderr, "computing the mean power spectrum of a random DNA sequence...\n");
    double dnaP[fMax+1];
    dnae_compute_random_dna_spectrum(N, dnaP);
    dnae_plot_spectra(dnaP, NULL, NULL, fMax, 0, dnaMax, FALSE, hSz, vSz, fSz, o->outName, 0, "-db");
    
    fprintf(stderr, "computing the power spectra of the two filters ...\n");
    wt_table_print(stderr, "initial", o->w0.ne, o->w0.e, 2);
    double w0P[fMax+1];
    dnae_compute_filter_spectrum(o->w0.ne, o->w0.e, N, w0P);
    dnae_plot_spectra(w0P, NULL, NULL, fMax, 0, wMax, FALSE, hSz, vSz, fSz, o->outName, 0, "-f0");

    wt_table_print(stderr, "incremental", o->w1.ne, o->w1.e, 2);
    double w1P[fMax+1];
    dnae_compute_filter_spectrum(o->w1.ne, o->w1.e, N, w1P);
    dnae_plot_spectra(w1P, NULL, NULL, fMax, 0, wMax, FALSE, hSz, vSz, fSz, o->outName, 0, "-f1");
    
    if (dnae_debug_fft) 
      { fprintf(stderr, "normalized power spectra of input data and filter kernels:\n");
        dnae_dump_spectra(fMax, "dnaP", dnaP, "w0P ", w0P, "w1P ", w1P);
      }

    fprintf(stderr, "running Fourier analysis up to level %d ...\n", o->maxLevel);
    
    /* !!! Should account for quantization error too !!! */
    int32_t f;
    
    /* Spectrum of filtered DNA signal, without any sampling: */
    double *fP = rn_alloc(fMax+1);
    for (f = 0; f < fMax; f++) { fP[f] = dnaP[f]; }
    
    /* Spectrum of non-aliased part of filtered and downsampled DNA signal: */
    double *vP = rn_alloc(fMax+1);
    for (f = 0; f < fMax; f++) { vP[f] = dnaP[f]; }
    
    /* Spectrum of aliased part of filtered and downsampled DNA signal: */
    double *aP = rn_alloc(fMax+1); 
    for (f = 0; f < fMax; f++) { aP[f] = 0; }
    
    /* Current filtering level: */
    int32_t level = 0;
    /* Current max freq in {vP,aP}: */
    int32_t gMax = fMax;
    while (TRUE)
      { /* Sampling step relative to base level: */
        int32_t step = (1 << level);
        fprintf(stderr, "step = %d gMax = %d\n", step, gMax);
        assert(gMax == N/step/2);
        
        /* Debugging: */
        if (dnae_debug_fft)
          { fprintf(stderr, "spectrum of filtered signal without sampling:\n");
            dnae_dump_spectra(fMax, "fP", fP, NULL, NULL, NULL, NULL);
            fprintf(stderr, "spectra of non-aliased and aliased parts of filtered sampled signals:\n");
            dnae_dump_spectra(gMax, "vP", vP, "aP", aP,   NULL, NULL);
          }
          
        /* Compute expected cutoff freq {hMax} for next level: */
        /* Plot spectrum of filtered, non-sampled signals: */
        /* Draw vertical line at the cutoff freq {gMax} of this level. */
        dnae_plot_spectra(fP, NULL, NULL, fMax, gMax, wMax, TRUE, hSz, vSz, fSz, o->outName, level, "-fo");
        
        /* Plot spectra of non-aliased and aliased filtered and sampled signals: */
        /* Draw vertical line at the cutoff freq {hMax} of the next level. */
        dnae_plot_spectra(vP, aP,   NULL, gMax, 0,    wMax, TRUE, hSz, vSz, fSz, o->outName, level, "-fs");

        /* Compute the current S,N powers and N/S ratio: */
        double vPt = 0; /* Total power in non-aliased freqencies. */
        double aPt = 0; /* Power in aliased frequencies. */
        int32_t g;
        for (g = 0; g < gMax; g++) { vPt += vP[g]; aPt += aP[g]; }
        fprintf
          ( stderr,
            "total power S+N = %8.5f  valid S = %8.5f  aliased N = %8.5f  N/(S+N) = %8.6f\n",
            vPt+aPt, vPt, aPt, aPt/(vPt+aPt)
          );

        /* Can we do another stage of filtering? */
        if (level >= o->maxLevel) { break; }
        if (gMax % 2 != 0) { fprintf(stderr, "!! odd sample length\n"); break; }
        if (o->w1.ne > N/step) { fprintf(stderr, "!! seq too short\n"); break; }

        /* Compute max frequency {hMax} after downsampling: */
        int32_t hMax = gMax/2; 
        assert(hMax == N/step/4);

        if (level == 0)
          { assert(gMax == fMax);
            for (g = 0; g <= gMax; g++) { fP[g] *= w0P[g]; vP[g] *= w0P[g]; aP[g] *= w0P[g]; }
          }
        else
          { /* Build effective kernel {wt[0..N-1]} of incremental filter in original domain: */
            double wt[N];
            int32_t t;
            for (t = 0; t < N; t++) { wt[t] = 0; }
            assert(o->w1.ne*step <= N);
            for (t = 0; t < o->w1.ne; t++) { wt[t*step] = o->w1.e[t]; }
            /* Compute power spectrum of incremental filter in original freq domain: */
            double wactP[fMax+1];
            dnae_compute_filter_spectrum(N, wt, N, wactP);
            /* Apply incremental filter to unsampled signal {fP}: */
            for (f = 0; f <= fMax; f++) { fP[f] *= wactP[f]; }
          
            /* Compute power spectrum of incremental filter in current freq domain: */
            double wincP[gMax+1];
            dnae_compute_filter_spectrum(o->w1.ne, o->w1.e, N/step, wincP);
            dnae_plot_spectra(wincP, NULL, NULL, gMax, hMax, wMax, TRUE, hSz, vSz, fSz, o->outName, level+1, "-fr");
            if (dnae_debug_fft)
              { fprintf(stderr, "power spectrum of incremental filter for level %d:\n", level+1);
                dnae_dump_spectra(fMax, "wincP", wincP, NULL, NULL, NULL, NULL);
              }

            /* Apply incremental filter to {vP,aP}: */
            for (g = 0; g <= gMax; g++) { vP[g] *= wincP[g]; aP[g] *= wincP[g]; }
          }


        /* Allocate new storage for valid and aliased spectra after downsampling: */
        double *vfP =  rn_alloc(hMax+1);
        double *afP =  rn_alloc(hMax+1);
        
        /* Simulate downsampling of {vP,aP} to {vfP,afP}: */
        int32_t h;
        for (h = 0; h <= hMax; h++) { vfP[h] = afP[h] = 0; }
        for (g = 0; g <= gMax; g++)
          { /* Compute the frequency {h} that {g} becomes when 2-sampled: */
            h = 2*g;
            if (h > gMax) h = h - 2*gMax;
            if (h < 0) { h = -h; }
            assert(h % 2 == 0);
            h /= 2;
            assert((h >= 0) && (h <= hMax));
            /* Add power {vP[g]}  to {vfP[g]} or {afP[g]} as appropriate: */
            if (g == h) 
              { /* Freq {g} gets preserved: */      vfP[h] += vP[g]; }
            else
              { /* Freq {g} gets aliased to {h}: */ afP[h] += vP[g]; }
            /* Add the power that was aliased in previous levels: */
            afP[h] += aP[g];
          }
          
        /* Prepare for next iteration: */
        free(vP); vP = vfP;
        free(aP); aP = afP;
        gMax = hMax;
        level++;
      }

    /* ...and we are done: */
    fprintf(stderr, "done.\n");
    return 0;
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

void dnae_compute_random_dna_spectrum(int32_t N, double dnaP[])
  { /* Get and check some parameters: */
    demand(N % 2 == 0, "N must be even");
    int32_t fMax = N/2; /* Maximum frequency n power spectrum. */
    int32_t nc = dnae_CHANNELS; /* Number of channels in datums. */
    assert(nc == 3); /* The analysis below is specific for {nc=3}. */
    double uP = 1.0; /* Expected power per channel and per freq {0..N-1} */
    /* In each channel, the samples are {-1} or {+1} with equal prob. */
    /* The expected value of the mean squared is {1/N}: */
    dnaP[0] = nc*uP;
    /* For frequency {fMax}, the expected power is {1/N} too: */
    dnaP[fMax] = nc*uP;
    /* For other freqs {1..fMax-1}, the expected power is {2/N}  too: */
    int32_t f;
    for (f = 1; f < fMax; f++) { dnaP[f] = 6*uP; }
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
  ( double *xP,
    double *yP,
    double *zP,
    int32_t fMax,
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
    asprintf(&levTag, "-%02d%s", level, tag);
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
    
    /* Set client reference window, with some skosh: */
    msm_ps_tools_set_graph_ref_window(dp, xMin, xMax, yMin, yMax);
    
    /* Shrink device window to leave space for tick labels: */
    msm_ps_tools_shrink_epswr_ref_window(dp, mrgL, mrgR, mrgB, mrgT);
    
    /* Decide whether to plot dots at individual samples: */
    double hStep = msm_ps_tools_map_x(dp, 1) - msm_ps_tools_map_x(dp, 0);
    bool_t show_dots = (hStep >= 1.5);  /* If spaced at least 1.5 mm */

    /* Draw axes, ticks, labels, etc: */
    char *font = "Times-Roman";
    double tickSize = fontSize/10.0;  /* In mm. */
    epswr_set_label_font(eps, font, fontSize);
    /* msm_ps_tools_draw_ref_axis(dp, HOR, 0.5,0.5,0.5); */
    /* msm_ps_tools_draw_ref_axis(dp, VER, 0.5,0.5,0.5); */
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
    
    argparser_get_keyword(pp, "-maxLevel");
    o->maxLevel = (int32_t)argparser_get_next_int(pp, 0, 20);

    argparser_get_keyword(pp, "-numPoints");
    o->numPoints = (int32_t)argparser_get_next_int(pp, 0, (1 << 18));

    argparser_get_keyword(pp, "-initFilter");
    o->w0 = wt_table_args_parse(pp, TRUE);

    argparser_get_keyword(pp, "-incrFilter");
    o->w1 = wt_table_args_parse(pp, TRUE);

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
    
    /* o->seqName = argparser_get_next(pp); */
    
    o->outName = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }

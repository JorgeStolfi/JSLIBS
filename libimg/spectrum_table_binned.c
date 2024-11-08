/* See spectrum_table_binned.h */
/* Last edited on 2024-11-06 10:11:03 by stolfi */ 

#define _GNU_SOURCE
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <jsmath.h>
#include <jsroot.h>
#include <affirm.h>
#include <bool.h>
#include <vec.h>
#include <float_image.h>

#include <spectrum_table_binned.h>

/* INTERNAL PROTOTYPES */

void spectrum_table_binned_freq_range(int32_t fr[], int32_t sz[], double *frLo, double *frHi);
  /* Given the frequency vector (in waves per image) and the image
    dimensions (in pixels) of a Hartley transform component, computes
    an interval {*frLo,*frHi} for its absolute frequency (in waves per
    pixel), assuming that its power is not concentrated at that
    frequency vector but spread over a small region of frequency space
    centered there.
    
    Specifically, the region is assumed to be an axis-aligned ellipse
    whose diameter along axis {i} is {2/sz[i]}, centered at
    coordinate {fr[i]/sz[i]} (both in waves per pixel). The returned interval
    is the orthogonal projection of that ellipse on the ray that 
    starts at {0,0} and goes through the center of that ellipse.
    
    As a special case, if {fr} is the null vector {(0,0)}, the 
    returned interval is the average for a ray with random direction. */

double spectrum_table_binned_compute_rel_terms_from_freq(double f);
  /* Computes the fraction of the Hartley terms that have frequency
    less than {f} waves per pixel, in the limit of large images.
  
    The frequency of a Hartley term is {hypot(fX,fY)} where {fX} and
    {fY} are unformly distributed over {[-0.5 _ +0.5]}. Therefore, the
    fraction of those terms with frequency less than {f} is the area
    of the intersection of {S} and {D(f)}, where {S} is the square of
    side 1 and {D(f)} is the disk of radius {f} concentric with
    {S}. */

double spectrum_table_binned_compute_freq_from_rel_terms(double a);
  /* The functional inverse of {spectrum_table_compute_rel_terms_from_freq}. */

double spectrum_table_binned_splat_weight(double frLo, double frHi, double fmin, double fmax);
  /* Computes the fraction of some term that should be splatted into a
    bucket whose frequency range is {[fmin _ fmax]}. Assumes that the
    term is spread over the frequency range {[fLo_fHi]} with a
    bell-shaped distribution. */

/* IMPLEMENTATIONS */

vec_typeimpl(spectrum_table_binned_t,spectrum_table_binned,spectrum_table_binned_entry_t); 
  
spectrum_table_binned_t spectrum_table_binned_make(int32_t nRanges)
  {
    demand(nRanges > 0, "invalid nRanges");
    int32_t i;
   /* Allocate binned table and define its fmid ranges: */
    spectrum_table_binned_t tb = spectrum_table_binned_new(nRanges);
    double fmin = 0; /* Lower limit of next band. */
    for (i = 0; i < nRanges; i++) 
      { spectrum_table_binned_entry_t *tbi = &(tb.e[i]);
        tbi->fmin = fmin;
        if (i == nRanges-1)
          { tbi->fmax = M_SQRT1_2; }
        else
          { /* Compute {tbi->fmax} so that the band contains {~1/nRanges} of all terms: */
            double rel_area_fmax = ((double)i + 1)/((double)nRanges);
            tbi->fmax = spectrum_table_binned_compute_freq_from_rel_terms(rel_area_fmax);
            /* Paranoia: */
            double r = spectrum_table_binned_compute_rel_terms_from_freq(tbi->fmax);
            demand(fabs(r - rel_area_fmax) < 1.0e-5, "bug in range comp");
          }
        assert(tbi->fmax > tbi->fmin);
        assert(tbi->fmax <= M_SQRT1_2);
        /* Compute {tbi->fmid} as the approximate median of the frequencies in band: */
        double rel_area_fmid = ((double)i + 0.5)/((double)nRanges);
        tbi->fmid = spectrum_table_binned_compute_freq_from_rel_terms(rel_area_fmid);
        assert(tbi->fmin <= tbi->fmid);
        assert(tbi->fmid <= tbi->fmax);
        /* Clear the accumulators: */
        tbi->nTerms = 0;
        tbi->power = 0;
        /* Prepare for next entry: */
        fmin = tbi->fmax;
      }
    return tb;
  }

void spectrum_table_binned_add_all
  ( float_image_t *P,  
    int32_t c,
    bool_t center,
    spectrum_table_binned_t *tb,
    bool_t verbose
  )
  {
    int32_t cols = (int32_t)P->sz[1];
    int32_t rows = (int32_t)P->sz[2];

    int32_t nRanges = tb->ne;
    demand(nRanges > 0, "invalid nRanges");
    
    if (verbose)
      { fprintf(stderr, "splatting spectrum onto %d bins\n", nRanges); }
 
    /* Total {nTerms} and {power} in channel {c} of {P}: */
    double P_nTerms = 0;
    double P_power = 0;
    
    /* General parameters of frequencies: */
    int32_t fd[2] = { cols, rows };        /* Denominator of natural freq vector. */
    int32_t fnMax[2] = { cols/2, rows/2 }; /* Max numerators of natural freq vector. */
    
    /* Splat the spectrum terms over the table {tb}, with smoothing. */
    for (int32_t ry = 0; ry < rows; ry++)
      for (int32_t rx = 0; rx < cols; rx++)
        { /* Adjust for centering: */
          int32_t fx = (center ? (rx + cols - cols/2) % cols : rx);
          int32_t fy = (center ? (ry + rows - rows/2) % rows : ry);
        
          /* Grab the spectrum entry's power: */
          double pwr = float_image_get_sample(P, c, fx, fy);

          /* Get its natural frequency vector: */
          int32_t fn[2]; /* Numerator of natural freq vector in waves/pixel. */
          fn[0] = (fx <= fnMax[0] ? fx : fx - cols);
          fn[1] = (fy <= fnMax[1] ? fy : fy - rows);
          
          /* Splat it onto {tb}: */
          spectrum_table_binned_add_term(tb, fn, fd, 1.0, pwr, verbose);
          
          /* Accumulate it into {P_nTerms,P_power}: */
          P_nTerms += 1.0;
          P_power += pwr;
        }
      
    if (verbose) 
      { /* Print total terms and power that were splatted onto {tb}: */
        fprintf(stderr, "added terms = %8.2f ", P_nTerms);
        fprintf(stderr, " power = %14.8f", P_power);
        fprintf(stderr, "\n");
      }
  }
  
void spectrum_table_binned_add_term
  ( spectrum_table_binned_t *tb,
    int32_t fn[], 
    int32_t fd[],
    double nTerms,
    double power,
    bool_t verbose
  )
  {
    /* Compute the signed range of frequencies spanned by the smoothed entry: */
    double frLo, frHi;
    spectrum_table_binned_freq_range(fn, fd, &frLo, &frHi);
    if (verbose) { fprintf(stderr, "  splatting  fr = [ %18.10f %18.10f ]\n", frLo, frHi); }
    /* The range {[frLo _ frHi]} must be folded over the range {H=[0.0 _ 0.707]}. */
    /* Compute the init reduced frequency {fr} and the direction {dir} for {frLo}: */
    double fr = frLo; /* Section {[frLo _ fr]} of splat done, section {[fr _ frHi]} remaining. */
    double fr_folded = fr - M_SQRT2*floor(fr/M_SQRT2);
    int32_t dir = +1;
    assert((fr_folded >= -1.0e-8) && (fr_folded < M_SQRT2+1.0e-8));
    if (fr_folded >= M_SQRT1_2) { fr_folded = M_SQRT2 - fr_folded; dir = -dir; }
    /* Fudge roundoff errors: */
    if (fr_folded < 0) { fr_folded = 0; }
    if (fr_folded > M_SQRT1_2) { fr_folded = M_SQRT1_2; }

    /* Locate the entry in table {tb} that contains {fr_folded}: */
    int32_t i = spectrum_table_binned_locate_entry(tb, fr_folded);
    while (fr < frHi)
      { if (verbose) { fprintf(stderr, "    fr = %18.10f  fr_folded = %18.10f", fr, fr_folded); }
        /* At this point, entry {tb.e[i]} contains {fr_folded}. */
        /* Also {fr_folded} is {fr} folded over {H}. */
        /* Also {dir} is the sign of the derivative of {fr_folded} wrt {fr}. */
        spectrum_table_binned_entry_t *tbi = &(tb->e[i]);
        assert(tbi->fmin <= fr_folded);
        assert(fr_folded <= tbi->fmax);
        /* Map {tbi}'s range from folded to unfolded frequency axis: */
        double df = fr - fr_folded;
        double sf = fr + fr_folded;
        double bin_fmin = (dir > 0 ? df + tbi->fmin : sf - tbi->fmax);
        double bin_fmax = (dir > 0 ? df + tbi->fmax : sf - tbi->fmin);
        if (verbose) { fprintf(stderr, "  bin = [ %18.10f  %18.10f ]\n", bin_fmin, bin_fmax); }
        /* Make sure that the portion already splatted will be excluded: */
        if (fr != frLo) { fr = bin_fmin; }
        /* Compute fraction of entry {txk} that falls inside the slice: */
        double w = spectrum_table_binned_splat_weight(frLo, frHi, bin_fmin, bin_fmax);
        /* Splat the entry: */
        tbi->nTerms += w * nTerms;
        tbi->power += w * power;
        /* Advance to the next bin (folding over the ends of {H}): */
        double fr_step;
        if (dir > 0)
          { fr_step = tbi->fmax - fr_folded;
            fr_folded = tbi->fmax;
            if (i == tb->ne - 1) { dir = -1; } else { i++; }
          }
        else
          { fr_step = fr_folded - tbi->fmin;
            fr_folded = tbi->fmin; 
            if (i == 0) { dir = +1; } else { i--; }
          }
        fr += fr_step; 
      }
    if (verbose) { fprintf(stderr, "  done splatting\n"); }
  }

void spectrum_table_binned_freq_range(int32_t fr[], int32_t sz[], double *frLo, double *frHi)
  {
    /* Compute the natural frequencies and splay widths on each axis (waves/pixel): */ 
    double fr_wpp[2]; /* Natural frequency in waves per pixel. */
    double rd_wpp[2]; /* Radius of fuzz ellipse along each axis. */
    int32_t i;
    for (i = 0; i < 2; i++)
      { /* Reduce the integer frequency modulo the image size: */
        int32_t szi = sz[i];
        int32_t frMax = szi/2;
        int32_t frNat = ((fr[i] % szi) + szi) % szi;
        if (frNat > frMax) { frNat = frNat - szi; }
        assert((-frMax <= frNat) && (frNat <= +frMax)); 
        /* Compute the axial natural frequency and frequency fuzz in waves per pixel: */
        fr_wpp[i] = ((double)frNat)/((double)szi);
        rd_wpp[i] = 1.0/((double)szi);
      }
    /* Compute the absolute natural frequency {fr_abs_wpp}: */
    double fr_abs_wpp = hypot(fr_wpp[0], fr_wpp[1]);
    /* Compute the splat width along the direction of {fr_wpp}: */
    double ux, uy;
    if (fr_abs_wpp == 0)
      { ux = uy = M_SQRT1_2; }
    else
      { ux = fabs(fr_wpp[0])/fr_abs_wpp;
        uy = fabs(fr_wpp[1])/fr_abs_wpp;
      }
    double rd_abs_wpp = hypot(rd_wpp[0]*ux, rd_wpp[1]*uy);
    /* Compute {*frLo,*frHi}: */
    (*frLo) = fr_abs_wpp - rd_abs_wpp;
    (*frHi) = fr_abs_wpp + rd_abs_wpp;
    /* fprintf(stderr, "  fr = ( %d/%d %d/%d )", fr[0], sz[0], fr[1], sz[1]); */
    /* fprintf(stderr, " = ( %10.8f %10.8f )", fr_wpp[0], fr_wpp[1]); */
    /* fprintf(stderr, " = %10.8f", fr_abs_wpp); */
    /* fprintf(stderr, " -> [ %10.8f  %10.8f ]", (*frLo), (*frHi)); */
    /* fprintf(stderr, "\n"); */
  }

int32_t spectrum_table_binned_locate_entry(spectrum_table_binned_t *tb, double f)
  {
    int32_t ntb = tb->ne;
    if (f <= tb->e[0].fmin) 
      { return 0; }
    else if (f >= tb->e[ntb-1].fmax)
      { return ntb-1; }
    else
      { int32_t iLo = 0;
        int32_t iHi = ntb-1;
        while (TRUE)
          { demand(iLo <= iHi, "table has gaps or is out of order");
            int32_t iMd = (iLo + iHi)/2;
            assert((iLo <= iMd) && (iMd <= iHi));
            if (f < tb->e[iMd].fmin)
              { iHi = iMd - 1; }
            else if (f > tb->e[iMd].fmax)
              { iLo = iMd + 1; }
            else
              { return iMd; }
          }
      }        
  }

double spectrum_table_binned_splat_weight(double frLo, double frHi, double fmin, double fmax)
  {
    if (fmin < frLo) { fmin = frLo; }
    if (fmax > frHi) { fmax = frHi; }
    
    if (fmin >= fmax) { return 0; }
    
    auto double sigmoid(double fr);
      /* A C2 sigmoid that goes from 0 to 1 over the interval {[frLo _ frHi]}. */
  
    double sigmoid(double fr)
      { double fx = (fr - frLo)/(frHi - frLo);
        return fx - sin(2*M_PI*fx)/(2*M_PI);
      }
      
    return sigmoid(fmax) - sigmoid(fmin);
  }

double spectrum_table_binned_compute_freq_from_rel_terms(double a)
  {
    if (a < M_PI/4)
      { return sqrt(a/M_PI); }
    else if (a >= 1.0)
      { return M_SQRT1_2; }
    else
      { 
        auto double G(double aa);
        /* An approximate inverse of {spectrum_table_binned_compute_rel_terms_from_freq}
           for {aa} in {[M_PI/4 _ 1.0]}. */
        
        double G(double aa) 
          { if (aa > 1) 
              { return 0; }
            else
              { return - sqrt((1 - aa)/(1 - M_PI/4)); }
          }
          
        double g = G(a); /* Desired value of {G(a)}. */

        auto double error(double freq); 
          /* A monotonic error measure between {freq} and the desired frequency: */
        
        double error(double freq)
          { double area = spectrum_table_binned_compute_rel_terms_from_freq(freq);
            return G(area) - g;
          }
      
        /* fprintf(stderr, "  looking for a = %18.10f g = %18.10f\n", a, g); */
        /* The current interval of {f} is {[fLo,fHi]}: */
        double frLo = 0.5;
        double frHi = M_SQRT1_2;
        return jsroot_bisect_secant(frLo, frHi, &error);
      }
  }
  
double spectrum_table_binned_compute_rel_terms_from_freq(double f)
  {
    double f2 = f*f;
    if (f < 0.5)
      { return M_PI * f2; }
    else if (f < M_SQRT1_2)
      { /* Area of Maltese cross: */
        double area_cross = 2*sqrt(f2 - 0.25); 
        /* Area of circular sectors betwen arms of Maltese cross: */
        double area_sects = f2*(M_PI - 4*acos(1/(2*f)));
        /* Total area: */
        return area_cross + area_sects;
      }
    else
      { /* Area of square: */
        return 1.0;
      }
  }


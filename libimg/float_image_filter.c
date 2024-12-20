/* See {float_image_filter.h}. */
/* Last edited on 2024-12-04 23:27:10 by stolfi */

#include <math.h>
#include <assert.h>
#include <values.h>
#include <string.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_filter.h>

/* INTERNAL PROTOTYPES */

void float_image_filter_dump_tables(FILE *wr, double s_lo[], double s_hi[], int32_t n);
  /* Writes the {s_lo[0..n-1]} and {s_hi[0..n-1]},
    assumed to be the weights of the low-freq and high-freq 
    low-pass filters for Filter frequencies {0..n-1}, 
    and the band-pass weights {s_hi[0..n-1] - s_lo[0..n-1]}. */

/* IMPLEMENTATIONS */

void float_image_filter_gaussian_band
  ( float_image_t *F, 
    r2_t *wMin, 
    r2_t *wMax, 
    bool_t complement,
    bool_t verbose
  )
  { int32_t chns = (int32_t)F->sz[0];
    int32_t cols = (int32_t)F->sz[1];
    int32_t rows = (int32_t)F->sz[2];
    
    /* Make the X and Y gaussian weight tables. */
    assert((0 <= wMin->c[0]) && (wMin->c[0] <= wMax->c[0]));
    assert((0 <= wMin->c[1]) && (wMin->c[1] <= wMax->c[1]));
    
    /* Weight tables that preserve terms with wavelength {> wMax}: */
    double *sx_fMin = float_image_filter_gaussian_freq_weights(cols, wMax->c[0]);
    double *sy_fMin = float_image_filter_gaussian_freq_weights(rows, wMax->c[1]);

    /* Weight tables that preserve terms with wavelength {> wMin}: */
    double *sx_fMax = float_image_filter_gaussian_freq_weights(cols, wMin->c[0]);
    double *sy_fMax = float_image_filter_gaussian_freq_weights(rows, wMin->c[1]);
    
    if (verbose)
      { fprintf(stderr, "# --- begin xweights.txt --------------------------------\n");
        float_image_filter_dump_tables(stderr, sx_fMin, sx_fMax, cols);
        fprintf(stderr, "# --- end xweights.txt --------------------------------\n");
        fprintf(stderr, "# --- begin yweights.txt --------------------------------\n");
        float_image_filter_dump_tables(stderr, sy_fMin, sy_fMax, rows);
        fprintf(stderr, "# --- end yweights.txt --------------------------------\n");
      }

    int32_t y, x, c;
    for (y = 0; y < rows; y++)
      { for (x = 0; x < cols; x++)
          { double s_fMin = sx_fMin[x]*sy_fMin[y];
            double s_fMax = sx_fMax[x]*sy_fMax[y];
            double s = s_fMax - s_fMin;
            if (complement) { s = 1 - s; }
            assert(s >= 0);
            for (c = 0; c < chns; c++)
              { float *smp = float_image_get_sample_address(F, c, x, y);
                (*smp) = (float)((*smp) * s);
              }
          }
      }
      
    free(sy_fMax);
    free(sx_fMax);
    free(sy_fMin);
    free(sx_fMin);
  }

double *float_image_filter_gaussian_freq_weights(int32_t n, double wRef)
  {
    double fRef = n/wRef;
    bool_t normSum = FALSE;
    bool_t folded = TRUE;
    double *w = gauss_table_make((uint32_t)n, 0.0, fRef, normSum, folded);
    /* If {wRef} is big but finite, returns the Dirac mask at freq {0,0}; */
    /* if {wRef} is {+INF}, returns a zero filter: */
    if (wRef == +INFINITY) { w[0] = 0.0; }
    return w;
  }
  
void float_image_filter_dump_tables(FILE *wr, double s_lo[], double s_hi[], int32_t n)
  {
    fprintf(wr, "#\n");
    fprintf(wr, "# %5s  %14s %14s  %14s\n", "F", "S_LO", "S_HI", "S"); 
    int32_t i;
    for (i = 0; i < n; i++)
      { double si = s_hi[i] - s_lo[i];
        fprintf(wr, "  %5d %14.11f %14.11f %14.11f\n", i, s_lo[i], s_hi[i], si);
      }
    fprintf(wr, "\n");
  }

/* See {float_image_hartley.h}. */
/* Last edited on 2017-01-02 21:31:45 by jstolfi */

#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <fftw3.h>
 
#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_hartley.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

void float_image_hartley_transform(float_image_t *A, float_image_t *T)
  { int chns = (int)A->sz[0];
    int cols = (int)A->sz[1];
    int rows = (int)A->sz[2];
    
    assert(chns == T->sz[0]);
    assert(cols == T->sz[1]);
    assert(rows == T->sz[2]);
    
    /* Allocate the work areas: */
    int N = (rows > cols ? rows : cols);
    double *in = (double*) fftw_malloc(sizeof(double) * N);
    double *out = (double*) fftw_malloc(sizeof(double) * N);
    int y, x, c;
    
    /* Do the row transforms: */
    fftw_plan px = fftw_plan_r2r_1d(cols, in, out, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    for (y = 0; y < rows; y++)
      { for (c = 0; c < chns; c++)
          { for (x = 0; x < cols; x++) { in[x] = float_image_get_sample(A, c, x, y); }
            fftw_execute(px);
            for (x = 0; x < cols; x++) { float_image_set_sample(T, c, x, y, (float)(out[x])); }
          }
      }
    fftw_destroy_plan(px);
    
    /* Do the column transforms: */
    fftw_plan py = fftw_plan_r2r_1d(rows, in, out, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    for (x = 0; x < cols; x++)
      { for (c = 0; c < chns; c++)
          { for (y = 0; y < rows; y++) { in[y] = float_image_get_sample(T, c, x, y); }
            fftw_execute(py);
            for (y = 0; y < rows; y++) { float_image_set_sample(T, c, x, y, (float)(out[y])); }
          }
      }
    fftw_destroy_plan(py);
    
    /* Combine elements to obtain pure sine-wave components: */
    int y0, x0;
    for (y0 = 1; y0 < (rows+1)/2; y0++)
      { int y1 = (rows - y0) % rows;
        for (x0 = 1; x0 < (cols+1)/2; x0++)
          { int x1 = (cols - x0) % cols;
            if ((x0 != x1) && (y0 != y1)) 
              { for (c = 0; c < chns; c++)
                  { float *smp00 = float_image_get_sample_address(T, c, x0, y0);
                    float *smp10 = float_image_get_sample_address(T, c, x1, y0);
                    float *smp01 = float_image_get_sample_address(T, c, x0, y1);
                    float *smp11 = float_image_get_sample_address(T, c, x1, y1);
                    double v00 = *smp00, v01 = *smp01, v10 = *smp10, v11 = *smp11;
                    *smp00 = (float)(+ v00 + v01 + v10 - v11)/2;
                    *smp10 = (float)(+ v00 - v01 + v10 + v11)/2;
                    *smp01 = (float)(+ v00 + v01 - v10 + v11)/2;
                    *smp11 = (float)(- v00 + v01 + v10 + v11)/2;
                  }
              }
          }
      }
    
    /* Scale elements to preserve sum of squares: */
    double s = sqrt(cols*rows);
    for (y = 0; y < rows; y++)
      { for (x = 0; x < cols; x++)
          { for (c = 0; c < chns; c++)
              { float *smp = float_image_get_sample_address(T, c, x, y);
                (*smp) = (float)((*smp)/s);
              }
          }
      }
        
    fftw_free(out);
    fftw_free(in);
  }

void float_image_hartley_wave(float_image_t *A, int fx, int fy, double amp)
  {
    int chns = (int)A->sz[0];
    int cols = (int)A->sz[1];
    int rows = (int)A->sz[2];
    
    /* Reduce the frequencies to the image domain: */
    fx = fx % cols; if (fx < 0) { fx += cols; }
    fy = fy % rows; if (fy < 0) { fy += rows; }
    assert((0 <= fx) && (fx < cols));
    assert((0 <= fy) && (fy < rows));
    
    int y, x, c; 
    for (y = 0; y < rows; y++)
      { for (x = 0; x < cols; x++)
          { for (c = 0; c < chns; c++)
              { float *smp = float_image_get_sample_address(A, c, x, y);
                double xCount = fx*((double)x)/((double)cols);
                double yCount = fy*((double)y)/((double)rows);
                (*smp) = (float)(amp*sin(M_PI*(2*xCount + 2*yCount + 0.25)));
              }
          }
      }
  }


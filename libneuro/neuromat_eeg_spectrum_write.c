/* See {neuromat_eeg_spectrum_write.h}. */
/* Last edited on 2023-12-14 08:33:57 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_eeg_spectrum_write.h>

void neuromat_eeg_spectrum_write
  ( FILE *wr, 
    int32_t nt, 
    int32_t ne, 
    int32_t kfmax, 
    double fsmp, 
    double **pwr
  )
  {
    demand(2*kfmax <= nt, "{kfmax} too high");
    
    for (int32_t kf = 0; kf <= kfmax; kf++) 
      { double f = kf*fsmp/nt;
        double flo = (kf == 0 ? 0.0 : kf - 0.5)*fsmp/nt;
        double fhi = (2*kf == nt ? (double)kf : kf + 0.5)*fsmp/nt;
        fprintf(wr, "%8d  %12.7f  %12.7f %12.7f ", kf, f, flo, fhi);
        for (int32_t ie = 0; ie < ne; ie++) { fprintf(wr, " %12.5e", pwr[ie][kf]); }
        fprintf(wr, "\n");
      }
    fflush(wr);
  }


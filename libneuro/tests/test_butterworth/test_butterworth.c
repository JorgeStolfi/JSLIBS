/* Tests the Butterworth filter formulas. */
/* Last edited on 2013-06-14 11:08:00 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#include <neuromat_filter.h>


int main(int argc, char **argv)
{

  int nf = 30000;     /* Number of Hartley/Fourier terms. */
  double fsmp = 300;  /* Assumed sampling frequency. */
  int npmax = 9;
  double fc[npmax+1];  /* Ref frequency for each order. */
  double fG2[npmax+1]; /* Freq where the order-{k} filter has gain 0.5. */
  double fGm[npmax+1]; /* Freq where the order-{k} filter has gain 0.001. */
    
  int kfmax = nf/2;
  int np, kf;
  for (np = 0; np <= npmax; np++) { fG2[np] = fGm[np] = NAN; fc[np] = 1.0; }
  for (kf = 1; kf < kfmax; kf++)
    { double f = (kf*fsmp)/nf;
      fprintf(stdout, "%10.4f", f);
      for (np = 1; np <= npmax; np++)
        { double G = sqrt(neuromat_filter_lowpass_butterworth(f, fc[np], np));
          fprintf(stdout, " %10.7f", G);
          if ((G < M_SQRT1_2) && (isnan(fG2[np]))) { fG2[np] = f; }
          if ((G < 0.001) && (isnan(fGm[np]))) { fGm[np] = f; }
        }
      fprintf(stdout, "\n");
    }
    
  fprintf(stderr, " %3s %10s %10s\n", "np", "fc/f2", "fc/fm"); 
  for(np = 1; np <= npmax; np++)
    { fprintf(stderr, " %3d %10.6f %10.6f\n", np, fc[np]/fG2[np], fc[np]/fGm[np]); }

  return 0;
}

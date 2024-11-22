/* See {tosl.h} */
/* Last edited on 2024-11-20 05:22:26 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include <haf.h>

#include <tosl.h>

double tosl_user_cpu_time_usec(void) 
  { struct timespec buf;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &buf);
    return (((double)buf.tv_sec)*1000000) + (((double)buf.tv_nsec)/1000);
  }

void tosl_compute_avg_dev(int32_t NT, double time[], double *avg_P, double *dev_P)
  { 
    double sum = 0;
    for (int32_t it = 0; it < NT; it++) { sum += time[it]; }
    double avg = sum/NT;
    double sum_dt2 = 0;
    for (int32_t it = 0; it < NT; it++) { double dt = time[it] - avg; sum_dt2 += dt*dt; }
    double dev = sqrt(sum_dt2/(NT-1));
    (*avg_P) = avg;
    (*dev_P) = dev;
  }

char *tosl_arc_id_to_string(tosl_arc_id_t ka)
  { char *res = NULL;
    if (ka == -1)
      { 
        char *res = jsprintf("**:*");
      }
    else
      { assert(ka >= 0);
        char *res = jsprintf("a%d:%d", ka/2, ka%2);
      }
    return res;
  }

void tosl_arc_id_print(FILE *wr, char *pref, tosl_arc_id_t ka, char *suff)
  { if (pref != NULL) { fputs(pref, wr); }
    char *xka = tosl_arc_id_to_string(ka);
    fputs(xka, wr);
    free(xka);
    if (suff != NULL) { fputs(suff, wr); }
  }

void tosl_tri_arc_id_print(FILE *wr, char *pref, tosl_arc_id_t ka0, tosl_arc_id_t ka1, tosl_arc_id_t ka2, char *suff)
  { if (pref != NULL) { fputs(pref, wr); }
    char *xka0 = tosl_arc_id_to_string(ka0);
    char *xka1 = tosl_arc_id_to_string(ka1);
    char *xka2 = tosl_arc_id_to_string(ka2);
    fprintf(wr, "%s - %s - %s", xka0, xka1, xka2);
    free(xka0); free(xka1); free(xka2);
    if (suff != NULL) { fputs(suff, wr); }
  }

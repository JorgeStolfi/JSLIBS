/* See {neuromat.h}. */
/* Last edited on 2013-05-30 23:58:23 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <fget.h>
#include <affirm.h>
#include <neuromat.h>

char **neuromat_e20_signal_names(void)
  { 
    int ne = 21;
    char **name = notnull(malloc(ne*sizeof(char*)), "no mem");
    name[ 0] = "F7";
    name[ 1] = "T3";
    name[ 2] = "T5";
    name[ 3] = "Fp1";
    name[ 4] = "F3";
    name[ 5] = "C3";
    name[ 6] = "P3";
    name[ 7] = "O1";
    name[ 8] = "F8";
    name[ 9] = "T4";
    name[10] = "T6";
    name[11] = "Fp2";
    name[12] = "F4";
    name[13] = "C4";
    name[14] = "P4";
    name[15] = "O2";
    name[16] = "Fz";
    name[17] = "Cz";
    name[18] = "Pz";
    name[19] = "Oz"; 
    name[20] = "TR"; /* Trigger(?) signal. */
    return name;
  }

double **neuromat_read_eeg_signals(FILE *rd, int nskip, int nt, int nc)
  { 
    double **val = notnull(malloc(nt*sizeof(double*)), "no mem");
    
    auto void skip_line(void);
      /* Consumes the remaining characters up to and including the EOL. 
        Bombs out if EOF before EOL. */
    
    int nl = 1; /* Internal line number (including comments and headers), from 1. */
    int nframes = 0; /* Data lines read (excluding comments and headers). */
    int nsaved = 0; /* Data lines stored (excluding skipped ones). */
    while (nsaved < nt) 
      { /* Read one more line: */
        int r = fgetc(rd);
        if (r == EOF)
          { fprintf(stderr, "** premature EOF at line %d\n", nl); exit(1); } 
        else if ((r == '#') || ((r >= 'A') && (r <= 'Z')) || ((r >= 'a') && (r <= 'z')))
          { /* Looks like a header line or comment: */
            if (nframes > 0) { fprintf(stderr, "** spurious header line at line %d\n", nl); exit(1); } 
            skip_line();
          }
        else if (nframes < nskip) 
          { skip_line(); nframes++; }
        else
          { ungetc(r, rd);
            double *valt = notnull(malloc(nc*sizeof(double)), "no mem");
            int ic;
            for (ic = 0; ic < nc; ic++) { valt[ic] = fget_double(rd); }
            (void)fget_test_char(rd, '\015');
            fget_eol(rd);
            val[nsaved] = valt;
            nsaved++;
            nframes++;
          }
        nl++;
      }
    return val;
    
    void skip_line(void)
      { 
        int r;
        do { r = fgetc(rd); } while ((r != EOF) && (r != '\n'));
        if (r == EOF) { fprintf(stderr, "** EOF found while skipping line %d\n", nl); exit(1); } 
      }
  }

void neuromat_write_eeg_signals(FILE *wr, int nt, int nc, double **val, int it_ini, int it_fin, int it_step)
  {
    demand((0 <= it_ini) && (it_ini <= it_fin) && (it_fin < nt), "invalid sample index range");
    demand(it_step > 0, "invalid sample index increment");
    int it, ic;
    int nw = 0;
    for (it = it_ini; it <= it_fin; it += it_step) 
      { for (ic = 0; ic < nc; ic++) { fprintf(wr, " %15.7e", val[it][ic]); }
        fprintf(wr, "\n");
        nw++;
      }
    fprintf(stderr, "wrote %d data frames\n", nw);
    fflush(wr);
  }

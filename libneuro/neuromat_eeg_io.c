/* See {neuromat_eeg_io.h}. */
/* Last edited on 2021-08-21 12:36:40 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <fget.h>
#include <affirm.h>
#include <jsstring.h>

#include <neuromat_eeg.h>

#include <neuromat_eeg_io.h>

int32_t neuromat_eeg_frame_read(FILE *rd, int32_t nc, double frm[], int32_t *nlP, int32_t *nfP)
  { int32_t nl = (nlP == NULL ? 0 : (*nlP)); /* Number of file lines already read/skipped (incl. comments and headers). */
    int32_t nf = (nfP == NULL ? 0 : (*nfP)); /* Number of data frames already read/skipped. */
    int32_t nfr = 0; /* Number of data lines read/skipped in this call. */
    
    auto void skip_line(void);
      /* Consumes the remaining characters up to and including the EOL. 
        Bombs out if EOF before EOL. */
    
    /* Loop until EOF or one data line: */
    while (TRUE)
      { /* Try to read one more line from the file: */
        int32_t r = fgetc(rd);
        if (r == EOF)
          { break; } 
        else
          { /* There is a line there: */
            nl++;
            if ((r == '#') || ((r >= 'A') && (r <= 'Z')) || ((r >= 'a') && (r <= 'z')))
              { /* Looks like a header line or comment, ignore it: */
                if (nf > 0) { fprintf(stderr, "!! spurious header/comment line at line %d, ignored\n", nl); } 
                skip_line();
              }
            else 
              { /* Looks like a data line: */
                nf++;
                nfr++;
                if (frm == NULL) 
                  { /* Just skip to EOL, inclusive: */
                    skip_line();
                  }
                else
                  { /* Parse the data values, store in {frm}: */
                    ungetc(r, rd);
                    for (int32_t ic = 0; ic < nc; ic++) { frm[ic] = fget_double(rd); }
                    (void)fget_skip_and_test_char(rd, '\015');
                    fget_eol(rd);
                  }
                break;
              }
          }
      }
      
    /* Update counters and return number of data lines read/skipped in this call: */
    (*nlP) = nl;
    (*nfP) = nf;
    return nfr;
    
    /* Internal implementations: */
    
    void skip_line(void)
      { int32_t r;
        do { r = fgetc(rd); } while ((r != EOF) && (r != '\n'));
        if (r == EOF) { fprintf(stderr, "** EOF found while skipping line %d\n", nl+1); exit(1); } 
      }
  }

double **neuromat_eeg_data_read(FILE *rd, int32_t nskip, int32_t nread, int32_t nc, int32_t *nlP, int32_t *ntP)
  { 
    int32_t nchunk = 3600*600;  /* One hour of samples at 600 Hz. */
    int32_t nalloc = (nread <= 0 ? nchunk : nread);  /* Trimmed or expanded later. */
    double **val = notnull(malloc(nalloc*sizeof(double*)), "no mem");
    
    int32_t nl = (nlP == NULL ? 0 : (*nlP)); /* Number of file lines already read/skipped (incl. comments and headers). */
    int32_t nf = 0; /* Number of data frames read/skiped (excluding comments and headers). */
    int32_t nt = 0; /* Number of data frames actually read and stored. */

    while ((nread <= 0) || (nt < nread))
      { 
        /* Try to read/skip one more data frame: */
        double *frm = (nf < nskip ? NULL : notnull(malloc(nc*sizeof(double)), "no mem"));
        int32_t nfr = neuromat_eeg_frame_read(rd, nc, frm, &nl, &nf);
        if (nfr == 0) 
          { /* Hit end of file: */
            if (frm != NULL) { free(frm); frm = NULL; }
            break;
          }
        else
          { /* Frame read/skip succeeded: */
            if (frm != NULL)
              { assert(nf > nskip);
                /* Save in frame list, expanding if needed: */
                if (nt >= nalloc)
                  { /* Allocation of {val} exhausted, expand: */
                    nalloc = nalloc + nchunk;
                    val = notnull(realloc(val, nalloc*sizeof(double*)), "no mem");
                  }
                val[nt] = frm;
                nt++;
              }
          }
      }
    assert(nf == nskip + nt); /* Paranoia. */
    if ((nread <= 0) || (nt < nread))
      { /* Stopped because of end-of-file. */
        fprintf(stderr, "reached end of input file after %d data frames\n", nf);
        if (nread > 0) { fprintf(stderr, "** expected at least %d frames\n", nskip+nread); exit(1); }
      }
    else
      { /* Stopped by {nread} limit. */
        assert(nt == nread);
      }
    fprintf(stderr, "skipped %d frames and parsed %d data frames\n", nf - nt, nt);
    if (nt < nalloc)
      { /* Trim excess storage in {val}: */
        val = realloc(val, nt*sizeof(double*));
        demand ((nt == 0) || (val != NULL), "no mem");
        nalloc = nt;
      }
    if (nlP != NULL) { (*nlP) = nl; }
    if (ntP != NULL) { (*ntP) = nt; }
    return val;
  }

void neuromat_eeg_data_write(FILE *wr, int32_t nt, int32_t nc, double **val, int32_t it_ini, int32_t it_fin, int32_t it_step)
  {
    demand((0 <= it_ini) && (it_ini <= it_fin) && (it_fin < nt), "invalid sample index range");
    demand(it_step > 0, "invalid sample index increment");
    int32_t nw = 0;
    for (int32_t it = it_ini; it <= it_fin; it += it_step) 
      { neuromat_eeg_frame_write(wr, nc, val[it]);
        nw++;
      }
    fprintf(stderr, "wrote %d data frames\n", nw);
    fflush(wr);
  }

void neuromat_eeg_frame_write(FILE *wr, int32_t nc, double val[])
  { for (int32_t ic = 0; ic < nc; ic++) { fprintf(wr, " %15.7e", val[ic]); }
    fprintf(wr, "\n");
  }

void neuromat_eeg_frame_print(FILE *wr, char *pre, int32_t nc, char **chnames, double val[], char *sep, char *suf)
  { if (pre != NULL) { fprintf(wr, "%s", pre); }
    for (int32_t i = 0; i < nc; i++)
      { if ((i > 0) && (sep != NULL)) { fprintf(wr, "%s", sep); } 
        if (chnames != NULL) { fprintf(wr, "%s = ", chnames[i]); }
        fprintf(wr, "%14.8e", val[i]);
      }
    if (suf != NULL) { fprintf(wr, "%s", suf); }
  }
 

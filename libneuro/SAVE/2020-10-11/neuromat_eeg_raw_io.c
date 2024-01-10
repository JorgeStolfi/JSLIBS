/* See {neuromat_eeg_raw_io.h}. */
/* Last edited on 2013-11-12 00:45:50 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>

#include <neuromat_eeg.h>

#include <neuromat_eeg_raw_io.h>

int neuromat_eeg_raw_frame_read(FILE *rd, int version, double unit, int nc, char **chnames, double val[])
  {
    int c = fgetc(rd);
    if (c == EOF) { return 0; }
    ungetc(c, rd);
    
    int i;
    for (i = 0; i < nc; i++)
      { val[i] = neuromat_eeg_raw_read_raw_sample(rd, version, unit, chnames[i]); }
    return 1;
  }

int16_t neuromat_eeg_raw_read_int16(FILE *rd, char *name, int32_t vmin, int32_t vmax)
  { 
    int c0 = fgetc(rd); demand(c0 != EOF, "unexpected eof: c0");
    int c1 = fgetc(rd); demand(c1 != EOF, "unexpected eof: c1");
    int16_t res = (int16_t)((c0 << 8) | c1);
    neuromat_eeg_raw_check_int_val(res, name, vmin, vmax);
    return res;
  }
  
int32_t neuromat_eeg_raw_read_int32(FILE *rd, char *name, int32_t vmin, int32_t vmax)
  {
    int c0 = fgetc(rd); demand(c0 != EOF, "unexpected eof: c0");
    int c1 = fgetc(rd); demand(c1 != EOF, "unexpected eof: c1");
    int c2 = fgetc(rd); demand(c2 != EOF, "unexpected eof: c2");
    int c3 = fgetc(rd); demand(c3 != EOF, "unexpected eof: c3");
    int32_t res = (((((c0 << 8) | c1) << 8) | c2) << 8) | c3;
    neuromat_eeg_raw_check_int_val(res, name, vmin, vmax);
    return res;
  }
  
float neuromat_eeg_raw_read_float(FILE *rd, char *name, float vmin, float vmax)
  {
    int c0 = fgetc(rd); demand(c0 != EOF, "unexpected eof: c0");
    int c1 = fgetc(rd); demand(c1 != EOF, "unexpected eof: c1");
    int c2 = fgetc(rd); demand(c2 != EOF, "unexpected eof: c2");
    int c3 = fgetc(rd); demand(c3 != EOF, "unexpected eof: c3");
    float res;
    uint32_t *p = (uint32_t *)(&res);
    (*p) = (((((c0 << 8) | c1) << 8) | c2) << 8) | c3;
    neuromat_eeg_raw_check_float_val(res, name, vmin, vmax);
    return res;
  }  

void neuromat_eeg_raw_check_int_val(int32_t res, char *name, int32_t vmin, int32_t vmax)
  {
    if ((res < vmin) || (res > vmax))
      { fprintf(stderr, "** field {%s} = %d should be in {%d .. %d}\n", name, res, vmin, vmax); 
        demand(FALSE, "** aborted");
      }
  }

void neuromat_eeg_raw_check_float_val(float res, char *name, float vmin, float vmax)
  {
    if ((res < vmin) || (res > vmax))
      { fprintf(stderr, "** field {%s} = %14.8e should be in {%14.8e .. %14.8e}\n", name, res, vmin, vmax); 
        demand(FALSE, "** aborted");
      }
  }

char *neuromat_eeg_raw_read_event_code(FILE *rd, char *name)
  {
    int nb = 4; /* Length of code in bytes. */
    char *res = notnull(malloc((nb+1)*sizeof(char)), "no mem");
    /* Read chars and store in {res}: */
    assert(nb == 4);
    int c0 = fgetc(rd); demand(c0 != EOF, "unexpected eof: c0");
    int c1 = fgetc(rd); demand(c1 != EOF, "unexpected eof: c1");
    int c2 = fgetc(rd); demand(c2 != EOF, "unexpected eof: c2");
    int c3 = fgetc(rd); demand(c3 != EOF, "unexpected eof: c3");
    res[0] = (char)c0;
    res[1] = (char)c1;
    res[2] = (char)c2;
    res[3] = (char)c3;
    res[4] = '\000';
    /* Remove all blanks: */
    int k = 0;
    int i;
    for (i = 0; i <= nb; i++) 
      { if (res[i] != ' ') 
          { res[k] = res[i]; k++; } 
      }
    return res;
  }

double neuromat_eeg_raw_read_raw_sample(FILE *rd, int version, double unit, char *name)
  {
    double res;
    if (version == 2)
      { int iv = neuromat_eeg_raw_read_int16(rd, name, -32767, +32767);
        res = unit*(double)iv;
      }
    else if (version == 4)
      { float fv = neuromat_eeg_raw_read_float(rd, name, -1000000.0f, +1000000.0f);;
        res = (double)fv;
      }
    else
      { demand(FALSE, "invalid raw file {version} code"); }
    demand(! isnan(res), "sample is NaN");
    return res;
  }

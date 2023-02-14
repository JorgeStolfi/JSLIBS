/* See {neuromat_eeg_source.h}. */
/* Last edited on 2023-02-12 07:52:33 by stolfi */
  
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include <fget.h>
#include <affirm.h>

#include <neuromat_eeg_header.h>

#include <neuromat_eeg_source.h>

neuromat_eeg_source_t *neuromat_eeg_source_new(void)
  {
    neuromat_eeg_source_t *ho = notnull(malloc(sizeof(neuromat_eeg_source_t)), "no mem");
    (*ho) = (neuromat_eeg_source_t)
      { .file = NULL,
        .nt = INT_MIN,
        .it_ini = INT_MIN,
        .it_fin = INT_MIN,
        .fsmp = NAN,
        .subject = INT_MIN,
        .run = INT_MIN       
      };
    return ho;
  }

void neuromat_eeg_source_free(neuromat_eeg_source_t *ho)
  {
    if (ho->file != NULL) { free(ho->file); }
    free(ho);
  }

void neuromat_eeg_source_write(FILE *wr, char *pref, neuromat_eeg_source_t *ho)
  {
    neuromat_eeg_header_write_field_string(wr, pref, "file", ho->file);
    neuromat_eeg_header_write_field_int(wr, pref, "nt", ho->nt, 1, INT_MAX-1);
    neuromat_eeg_header_write_field_int_range(wr, pref, "sample_range", ho->it_ini, ho->it_fin, INT_MIN+1, INT_MAX-1);
    neuromat_eeg_header_write_field_double(wr, pref, "fsmp", ho->fsmp, 0.01, 1.0e12);
    neuromat_eeg_header_write_field_int(wr, pref, "subject", ho->subject, 1, INT_MAX-1);
    neuromat_eeg_header_write_field_int(wr, pref, "run", ho->run, 1, INT_MAX-1);
    fflush(wr);
  }

void neuromat_eeg_source_read_field(FILE *rd, char *name, neuromat_eeg_source_t *ho)
  {
    if (strcmp(name, "file") == 0) 
      { ho->file = fget_string(rd); }
    else if (strcmp(name, "nt") == 0) 
      { ho->nt = fget_int32(rd); 
        demand(ho->nt > 0, "invalid original file length");
      }
    else  if (strcmp(name, "sample_range") == 0) 
      { ho->it_ini = fget_int32(rd); 
        ho->it_fin = fget_int32(rd); 
        demand(ho->it_fin >= ho->it_ini, "invalid sample index range");
      }
    else if (strcmp(name, "fsmp") == 0) 
      { ho->fsmp = fget_double(rd);
        demand((! isnan(ho->fsmp)) && (ho->fsmp > 0), "invalid sampling frequency");
      }
    else if (strcmp(name, "subject") == 0) 
      { ho->subject = fget_int32(rd); 
        demand(ho->subject > 0, "invalid subject number");
      }
    else if (strcmp(name, "run") == 0) 
      { ho->run = fget_int32(rd); 
        demand(ho->run > 0, "invalid run number");
      }
    else
      { demand(FALSE, "invalid source info field name"); }
  }

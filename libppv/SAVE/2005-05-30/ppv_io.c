/* See ppv_io.h */
/* Last edited on 2005-05-29 11:11:20 by stolfi */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#include <ppv_array.h>
#include <ppv_io.h>

#include <filefmt.h>

#include <stdlib.h>
#include <stdio.h>

/* INTERNAL PROTOTPES */

void ppv_write_pixels_plain ( FILE *wr, ppv_array_t *A );
void ppv_write_pixels_raw_small ( FILE *wr, ppv_array_t *A );
void ppv_write_pixels_raw_big ( FILE *wr, ppv_array_t *A );
void ppv_write_pixels_raw_bytes ( FILE *wr, ppv_array_t *A );

#define ppv_FILE_TYPE "ppv_array_t"
#define ppv_FILE_VERSION "2005-05-28"
    
void ppv_write_array(FILE *wr, ppv_array_t *A, bool_t plain)
{
  ppv_axis_t ax;
  
  /* Write header: */
  filefmt_write_header(wr, ppv_FILE_TYPE, ppv_FILE_VERSION);
  fprintf(wr, "nax = %d\n", ppv_NAX);
  fprintf(wr, "size = ");
  for (ax = 0; ax < ppv_NAX; ax++) { fprintf(wr, " %d", A->size[ax]); }
  fprintf(wr, "\n");
  fprintf(wr, "bps = %d", A->bps);
  fprintf(wr, "plain = %d", (plain ? 1 : 0));
  
  /* Now loop on samples, first index faster: */
  if (! plain)
    { ppv_write_pixels_plain(wr, A); }
  else if (A->bps < 8)
    { ppv_write_pixels_raw_small(wr, A); }
  else if (A->bps > 8)
    { ppv_write_pixels_raw_big(wr, A); }
  else /* A->bps == 8 */
    { ppv_write_pixels_raw_bytes(wr, A); }     
  fputc('\n', wr);

  /* Write header: */
  filefmt_write_footer(wr, ppv_FILE_TYPE);

  fflush(wr);
}

void ppv_write_pixels_plain ( FILE *wr, ppv_array_t *A )
  {
    ppv_index_t ix[ppv_NAX];
    for (ix[5] = 0; ix[5] < A->size[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < A->size[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < A->size[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < A->size[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < A->size[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < A->size[0]; ix[0]++)
                { ppv_pos_t pos = ppv_sample_pos(A, ix);
                  ppv_sample_t qv = ppv_get_sample(A->el, A->bps, A->bpw, pos);
                  fprintf(wr, " %u", qv);
                }                  
  }

void ppv_write_pixels_raw_small ( FILE *wr, ppv_array_t *A )
  {
    int spc = 8/A->bps; /* Samples per byte. */
    ppv_word_08_t buf = 0;
    ppv_nbits_t shift_ini = (spc - 1)*A->bps; /* Shift to apply to sample [0]. */
    ppv_nbits_t shift = shift_ini; /* Shift to apply to next sample. */
    ppv_index_t ix[ppv_NAX];
    for (ix[5] = 0; ix[5] < A->size[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < A->size[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < A->size[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < A->size[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < A->size[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < A->size[0]; ix[0]++)
                { ppv_pos_t pos = ppv_sample_pos(A, ix);
                  ppv_sample_t qv = ppv_get_sample(A->el, A->bps, A->bpw, pos);
                  buf = buf | (qv << shift);
                  if (shift == 0) 
                    { fputc((char)buf, wr); shift = shift_ini; buf = 0; }
                  else
                    { shift -= A->bps; }
                }
    if (shift < shift_ini) { fputc((char)buf, wr); }
  }

void ppv_write_pixels_raw_big ( FILE *wr, ppv_array_t *A )
  {
    ppv_nbits_t cps = (A->bps + 7)/8; /* Bytes per sample. */
    ppv_word_t maskw = 255; /* Mask to chop a byteful from the sample. */
    ppv_nbits_t shift_ini = (cps-1)*8; /* Shift to apply to sample to get first byte. */
    ppv_word_08_t buf;
    ppv_nbits_t shift = shift_ini; /* Shift to get next byte. */
    ppv_index_t ix[ppv_NAX];
    for (ix[5] = 0; ix[5] < A->size[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < A->size[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < A->size[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < A->size[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < A->size[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < A->size[0]; ix[0]++)
                { ppv_pos_t pos = ppv_sample_pos(A, ix);
                  ppv_sample_t qv = ppv_get_sample(A->el, A->bps, A->bpw, pos);
                  int k;
                  for (k = 1; k < cps; k++) 
                    { buf = (qv >> shift) & maskw; 
                      fputc((char)buf, wr); 
                      shift -= A->bps;
                    }
                } 
  }
                  
void ppv_write_pixels_raw_bytes ( FILE *wr, ppv_array_t *A )
  {
    ppv_index_t ix[ppv_NAX];
    for (ix[5] = 0; ix[5] < A->size[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < A->size[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < A->size[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < A->size[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < A->size[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < A->size[0]; ix[0]++)
                { ppv_pos_t pos = ppv_sample_pos(A, ix);
                  ppv_sample_t qv = ppv_get_sample(A->el, A->bps, A->bpw, pos);
                  fputc((char)qv, wr);
                }                  
  }

ppv_array_t ppv_read_array(FILE *rd, ppv_nbits_t bpw)
  {
    ppv_array_t A;
    
    /* Read header: */
    filefmt_read_header(rd, ppv_FILE_TYPE, ppv_FILE_VERSION);
    int nax, bps;
    
    nax = nget_int(rd, "nax"); fget_eol(rd);
    if (nax > ppv_NAX) { ppv_error("image file has too many indices"); }
    nget_name_eq(rd, "size");
    for (ax = 0; ax < nax; ax++) 
      { 
    
    
    
    /* Read footer: */
    filefmt_read_footer(rd, ppv_FILE_TYPE);
    return A;
  }


/* See {rdag_io.h}. */
/* Last edited on 2017-01-04 18:01:18 by stolfilocal */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <affirm.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <bool.h>

#include <rdag.h>
#include <rdag_def.h>
#include <rdag_io.h>

#define rdag_FILE_TYPE "rdag_t"
#define rdag_FILE_VERSION "2009-10-24"

void rdag_write(FILE *wr, rdag_t *D)
  {
    filefmt_write_header(wr, rdag_FILE_TYPE, rdag_FILE_VERSION);

    fprintf(wr, "nn = %u\n", D->nn);
    fprintf(wr, "ni = %u\n", D->ni);
    fprintf(wr, "no = %u\n", D->no);
    fprintf(wr, "max_node = %u\n", D->max_node);
    
    uint32_t k;
    for (k = 0; k < D->max_node; k++)
      { rdag_node_t s = k + 1;
        rdag_node_data_t dt;
        rdag_node_data_get(D, s, &dt);
        fprintf(wr, "%u %u %u %u %u\n", s, dt.f_link, dt.i_mark, dt.o_mark, dt.p_link);
      }
    filefmt_write_footer(wr, rdag_FILE_TYPE);
    fflush(wr);
  }

rdag_t *rdag_read(FILE *rd)
  {
    /* Read and check the file header: */
    filefmt_read_header(rd, rdag_FILE_TYPE, rdag_FILE_VERSION);
    
    /* Get the bit widths: */
    uint32_t nn = nget_uint32(rd, "nn", 10); fget_eol(rd);
    uint32_t ni = nget_uint32(rd, "ni", 10); fget_eol(rd);
    uint32_t no = nget_uint32(rd, "no", 10); fget_eol(rd);
    
    /* Create an empty automaton {D} with those bit widths: */
    rdag_t *D = rdag_new(nn, ni, no, 0);
    
    /* Read the max node number and try to expand {D} to that size: */
    uint32_t max_node = nget_uint32(rd, "max_node", 10); fget_eol(rd);
    rdag_expand(D, max_node);
    D->max_node = max_node;
    D->hash_valid = FALSE;
    
    /* Read the node data: */
    uint32_t k;
    for (k = 0; k < max_node; k++)
      { rdag_node_t s = fget_uint32(rd, 10);
        demand(s == k + 1, "nodes out of sequence");
        rdag_node_data_t dt;
        dt.f_link = fget_uint32(rd, 10);
        dt.i_mark = fget_uint32(rd, 10);
        dt.o_mark = fget_uint32(rd, 10);
        dt.p_link = fget_uint32(rd, 10);
        fget_eol(rd);
        rdag_node_data_set(D, s, &dt);
      }
    
    /* Read and check the file footer: */
    filefmt_read_footer(rd, rdag_FILE_TYPE);
    return D;
  }

int rdag_line_read (FILE *rd, uint32_t *nbuf, char buf[], uint32_t nmax)
  {
    int c;
    (*nbuf) = 0;

    /* Read first char of line (maybe space or end-of-line), stop if EOF: */
    c = fgetc(rd);
    if (c == EOF) { return -1; }
    
    uint64_t nread = 0; /* Local number of characters read. */
    while ((c != '\012') && (c != '\015'))
      { if (nread < (uint64_t)nmax) { buf[nread] = (char)c; (*nbuf)++;  }
        nread++;
        c = fgetc(rd);
        /* If the last line is not empty and lacks a end-of-line, provide one: */
        if (c == EOF) { c = '\012'; }
      }
    
    if (c == '\015')
      { /* Discard LF after CR: */
        c = fgetc(rd);
        if (c != '\012') 
          { /* Oops, CR only: */
            if (c != EOF) { ungetc(c, rd); }
          }
      }

    assert((*nbuf) < nmax);
    return (nread <= (uint64_t)nmax ? 0 : -2);
  }

void rdag_string_iso_latin_1_cleanup(uint32_t *nbuf, char buf[])
  { 
    uint32_t k = 0;  /* Input characters already processed. */
    uint32_t m = 0;  /* Byte count in cleaned string. */
    while (k < (*nbuf))
      { unsigned char c = buf[k]; k++;
        /* Convert whitespace: */
        if (((c >= '\011') && (c <= '\015')) || (c == (unsigned char)'\240')) 
          { c = ' '; }
        /* Store in output string: */
        if (c == ' ')
          { /* Store only if the preceding char was not space: */
            if ((m > 0) && (buf[m-1] != ' '))
              { buf[m] = c; m++; }
          } 
        else
          { buf[m] = c; m++; }
      }
    /* Delete the trailing NUL if any: */
    if ((m > 0) && (buf[m-1] == '\000')) { m--; }
    /* Delete the trailing space if any: */
    if ((m > 0) && (buf[m-1] == ' ')) { m--; }
    /* Update string length: */
    (*nbuf) = m;
  }

bool_t rdag_string_iso_latin_1_encode(uint32_t nbuf, char buf[], uint32_t *nstr, rdag_symbol_t str[], uint32_t nmax)
  { 
    (*nstr) = 0; 
    uint32_t k = 0;  /* Input characters already processed. */
    while (k < nbuf)
      { unsigned char c = buf[k]; k++;
        /* Check for bizarre spaces or non-printable chars: */
        if ((c <= '\037') || ((c >= '\177') && (c <= (unsigned char)'\240'))) 
          { return FALSE; }
        /* Check for soft hyphen: */
        if (c == (unsigned char)'\255') 
          { return FALSE; }
        /* Check max length: */
        if ((*nstr) >= nmax) { return FALSE; }
        /* Store symbol: */
        str[*nstr] = (rdag_symbol_t)c; (*nstr)++;
      }
    return TRUE;
  }

void rdag_string_iso_latin_1_write(FILE *wr, uint32_t nstr, rdag_symbol_t str[])
  {
    uint32_t k;
    for (k = 0; k < nstr;k++) { fputc((char)str[k], wr); }
  }

/* See {imq.h}  */
/* Last edited on 2024-12-05 01:09:43 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <imq.h>

int32_t imq_read_var_length_record(FILE *rd, int32_t nb, byte_t ibuf[])
  { 
    /* Read the record length, a 16-bit number, little-endian: */
    int32_t lo = getc(rd);
    if (lo == EOF) { return 0; }
    int32_t hi = getc(rd);
    demand(hi != EOF, "EOF in middle of record length");
    int32_t length = ((hi << 8) | lo);
    demand(length <= nb, "not enough space in {ibuf}");
    
    /* Try to read {length} bytes into {ibuf}: */
    int32_t nlen = (int32_t)fread(ibuf, 1, (uint32_t)length, rd);
    demand (nlen == length, "EOF in middle of record body");

    if ((length & 1) != 0)
      { /* Read the extra filler byte: */
        int32_t fill = getc(rd);
        demand(fill != EOF, "EOF when reading even-fill byte");
      }
    return (int32_t)length;
  }

void imq_read_header
  ( FILE *rd, 
    char *format_P, 
    int32_t *record_bytes_P, 
    int32_t *label_checksum_P,
    int32_t *NX_P, 
    int32_t *NY_P, 
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "reading input IMQ header...\n"); }

    int32_t nb_max = 2048;
    byte_t ibuf[nb_max];

    char format = '?';
    int64_t record_bytes = -1;
    int64_t label_checksum = -1;

    /* Read the PDS labels to the "END" or end-of-file: */
    while (TRUE)
      { int32_t length = imq_read_var_length_record(rd, nb_max, ibuf);
        if (length == 0) { break; }
        
        /* Isolate the keyword: */
        char *pbeg = (char *)ibuf; 
        char *pkey = pbeg; while (((pkey - pbeg) < length) && ((*pkey) == ' ')) { pkey++; } 
        char *pend = pkey; while (((pend - pbeg) < length) && ((*pend) != ' ')) { pend++; }
        int32_t nc = (int32_t)(pend-pkey);
        if (verbose) { fprintf(stderr, "read %.*s label record (%u bytes)\n", nc, pkey, length); }
        
        /* Check for useful keywords: */
        if ((nc == 8) && (strncmp(pkey, " CHECKSUM", 8) == 0))
          { /* Get the checksum and store in label_checksum: */
            demand(label_checksum == -1, "duplicate CHECKSUM field");
            int32_t nf = sscanf((char*)(ibuf+35), "%lu", &label_checksum);
            demand(nf == 1, "invalid CHECKSUM value");
          }
        else if ((nc == 12) && (strncmp(pkey, "RECORD_BYTES", 12) == 0))
          { demand(record_bytes == -1, "duplicate RECORD_BYTES field");
            int32_t nf = sscanf((char*)(ibuf+35), "%ld", &record_bytes);
            demand(nf == 1, "invalid RECORD_BYTES value");
            if ((record_bytes != 836) && (record_bytes != 1204))
              { fprintf(stderr, "imqtopgm: RECORD_BYTES is %ld, changed to %d\n", record_bytes, 1204);
                record_bytes = 1204;
              }
          }
        else if ((nc == 3) && (strncmp(pkey, "END", 3) == 0)) 
          { break; }
      } 

    if (verbose) { fprintf(stderr, "determining format and image size...\n"); }
    int32_t NX;
    int32_t NY;
    if (record_bytes == 836)
      { /* Voyager image: */
        format = 'Y'; NY =  800; NX = 800;
      }
    else
      { /* Viking image: */
        format = 'K'; NY = 1056; NX = 1204;
      }

    /* Return values to caller: */
    (*format_P) = format;
    (*record_bytes_P) = (int32_t)record_bytes;
    (*label_checksum_P) = (int32_t)label_checksum;
    
    (*NX_P) = NX;
    (*NY_P) = NY;
  }

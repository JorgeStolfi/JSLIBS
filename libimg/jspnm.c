/* See jspnm.h */
/* Last edited on 2023-11-26 07:10:18 by stolfi */

#define _GNU_SOURCE
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>

#include <affirm.h>
#include <jspnm.h>
#include <jsfile.h>
#include <fget.h>
#include <rn.h>

void pnm_choose_output_format
  ( uint16_t maxval, 
    int chns, 
    bool_t forceplain,
    pnm_format_t *formatP,
    bool_t *rawP,
    bool_t *bitsP
  )
  { pnm_format_t format;
    uint16_t maxraw = PNM_FILE_MAX_MAXVAL;
    if (maxval <= 0){ pnm_error("invalid maxval (%u)", maxval); }
    bool_t raw = (maxval <= maxraw) & (! forceplain);
    bool_t bits = (chns == 1) & (maxval == 1);
    if (bits)
      { format = (raw ? RPBM_FORMAT : PBM_FORMAT); }
    else if (chns == 1)
      { format = (raw ? RPGM_FORMAT : PGM_FORMAT); }
    else if (chns == 3)
      { format = (raw ? RPPM_FORMAT : PPM_FORMAT); }
    else
      { pnm_error( "bad channel count (%d)", chns); format = 0; }
    if (maxval > pnm_file_max_maxval(format))
      { pnm_error("invalid maxval (%u) for %c%c file", maxval, format/256, format%256); }
    *formatP = format;
    *rawP = raw;
    *bitsP = bits;
  }

void pnm_read_header
  ( FILE *rd, 
    int *colsP, 
    int *rowsP, 
    int *chnsP, 
    uint16_t *maxvalP, 
    bool_t *rawP, 
    bool_t *bitsP, 
    pnm_format_t *formatP
  )
  { /* Read and check magic number. */
    pnm_format_t format = pnm_read_format(rd);
    *formatP = format;
    
    /* Deduce number of channels {chns} and raw/plain flag from {format}: */
    int chns = 0;
    bool_t raw = FALSE;
    bool_t bits = FALSE;
    switch(format)
      { case PBM_FORMAT:  chns = 1; raw = FALSE; bits = TRUE;  break;
        case RPBM_FORMAT: chns = 1; raw = TRUE;  bits = TRUE;  break;
        case PGM_FORMAT:  chns = 1; raw = FALSE; bits = FALSE; break;
        case RPGM_FORMAT: chns = 1; raw = TRUE;  bits = FALSE; break; 
        case PPM_FORMAT:  chns = 3; raw = FALSE; bits = FALSE; break;
        case RPPM_FORMAT: chns = 3; raw = TRUE;  bits = FALSE; break;
        default: pnm_error( "bad format (%c%c) -- neither PBM/PGM/PPM", format/256, format%256 );
      }
    *chnsP = chns;
    *rawP = raw;
    *bitsP = bits;
    
    /* Skip comments and read width: */
    pnm_skip_whitespace_and_comments(rd, FALSE);
    *colsP = pnm_read_plain_int(rd);

    /* Skip comments and read height: */
    pnm_skip_whitespace_and_comments(rd, FALSE);
    *rowsP = pnm_read_plain_int(rd);
    
    uint16_t maxval;
    if ((format == PBM_FORMAT) || (format == RPBM_FORMAT))
      { /* PBM files do not have a {maxval} field in header: */
        maxval = 1;
      }
    else
      { /* Skip comments and read maxval: */
        uint16_t maxmaxval = pnm_file_max_maxval(format);
        pnm_skip_whitespace_and_comments(rd, FALSE);
        maxval = pnm_read_plain_sample(rd, maxmaxval);
      }
    *maxvalP = maxval;
    
    /* Skip ONE char that must be blank, newline, etc.: */
    int sp = fgetc(rd);
    if ((sp != '\000') && (sp != ' ') && ((sp < '\011') || (sp > '\015')))
      { pnm_error("EOF or bad character after maxval"); }
  }

void pnm_write_header(FILE *wr, int cols, int rows, uint16_t maxval, pnm_format_t format)
  { /* Write file's magic number: */
    pnm_write_format(wr, format);
    /* Write dimensions: */
    fprintf(wr, "\n%d %d\n", cols, rows);
    if ((format != PBM_FORMAT) && (format != RPBM_FORMAT))
      { /* Write the maxval field: */
        fprintf(wr, "%u\n", maxval);
      }
    fflush(wr);
  }

pnm_format_t pnm_read_format(FILE *rd)
  { int c1 = getc(rd); if (c1 == EOF) { pnm_error("empty image file"); }
    int c2 = getc(rd); if (c2 == EOF) { pnm_error("unexpected EOF in magic numbr"); }
    return (pnm_format_t)((c1 << 8) | c2);
  }

void pnm_write_format(FILE *wr, pnm_format_t format)
  { putc((char)((format/256)&255), wr);
    putc((char)(format&255), wr);
  }

void pnm_skip_whitespace_and_comments(FILE *rd, bool_t verbose)
  { int c = fgetc(rd);
    while ((c == '#') || (c == ' ') || (c == '\011') || (c == '\012') || (c == '\015')) 
      { if (c == '#')
          { do
              { if (verbose) { fprintf(stderr, "%c", c); }
                c = fgetc(rd);
                if (c == EOF) { return; }
              }
            while (c != '\n');
            if (verbose) { fprintf(stderr, "%c", c); }
          }
        c = fgetc(rd);
      }
    if (c == EOF) { return; }
    ungetc(c, rd);
  }

uint16_t pnm_file_max_maxval(pnm_format_t format)
  { uint16_t maxmax = 0;
    switch(format)
      { case PBM_FORMAT:  maxmax = 1; break;
        case RPBM_FORMAT: maxmax = 1; break;
        case PGM_FORMAT:  maxmax = PNM_FILE_MAX_MAXVAL; break;
        case RPGM_FORMAT: maxmax = PNM_FILE_MAX_MAXVAL; break; 
        case PPM_FORMAT:  maxmax = PNM_FILE_MAX_MAXVAL; break;
        case RPPM_FORMAT: maxmax = PNM_FILE_MAX_MAXVAL; break;
        default: 
          pnm_error( "bad format (%c%c) -- neither PBM/PGM/PPM", format/256, format%256 );
      }
    return maxmax;
  }

int pnm_read_plain_int(FILE *rd)
  { return fget_int32(rd); }

void pnm_write_plain_int(FILE *wr, int ival)
  { fprintf(wr, "%d", ival); }

uint16_t pnm_read_plain_sample(FILE *rd, uint16_t maxval)
  { /* Skip any blanks or newlines: */
    fget_skip_formatting_chars(rd);
    /* Read and check the sample value {ival}: */
    unsigned int x = fget_uint32(rd, 10);
    assert(x <= PNM_MAX_SAMPLE);
    uint16_t ival = (uint16_t)x;
    pnm_check_sample_range(&ival, maxval);
    return ival;
  }

int pnm_write_plain_sample(FILE *wr, uint16_t ival, uint16_t maxval)
  { pnm_check_sample_range(&ival, maxval);
    /* Extract digits of {ival} in {buf[0..n-1]}: */
    char buf[30];
    int n = 0;
    do { buf[n] = (char)('0' + (ival % 10)); ival /= 10; maxval /= 10; n++; } while (ival > 0); 
    /* Complete with blanks to same size as {maxval}: */
    while (maxval > 0) { buf[n] = ' '; maxval /= 10; n++; }; 
    /* Print digits in correct order: */
    int i = n;
    while (i > 0) { i--; putc(buf[i], wr); } 
    return n;
  }

uint16_t pnm_read_plain_bit(FILE *rd)
  { /* Skip any blanks or newlines: */
    fget_skip_formatting_chars(rd);
    /* Read and check the sample value: */
    int c1 = getc(rd); if (c1 == EOF) { pnm_error("unexpected EOF in sample"); }
    if (c1 == '0') 
      { return 0; }
    else if (c1 == '1')
      { return 1; }
    else
      { pnm_error("invalid sample in PBM file");
        return 0; 
      }
  }

void pnm_write_plain_bit(FILE *wr, uint16_t ival)
  { if (ival == 0)
      { fputc('0', wr); }
    else if (ival == 1)
      { fputc('1', wr); }
    else
      { pnm_error("invalid sample for PBM file"); }
  }

uint16_t pnm_read_raw_byte(FILE *rd, uint16_t maxval)
  { int c1 = getc(rd); if (c1 == EOF) { pnm_error("unexpected EOF in sample"); }
    uint16_t ival = (uint16_t)(c1 & 255);
    pnm_check_sample_range(&ival, maxval);
    return ival;
  }

void pnm_write_raw_byte(FILE *wr, uint16_t ival, uint16_t maxval)
  { pnm_check_sample_range(&ival, maxval);
    putc((char)(ival&255), wr);
  }

uint16_t pnm_read_raw_short(FILE *rd, uint16_t maxval)
  { int c1 = getc(rd); if (c1 == EOF) { pnm_error("unexpected EOF in sample"); }
    int c2 = getc(rd); if (c2 == EOF) { pnm_error("unexpected EOF in sample"); }
    uint16_t ival =  (uint16_t)((c1 & 255) * 256 + (c2 & 255));
    pnm_check_sample_range(&ival, maxval);
    return ival;
  }

void pnm_write_raw_short(FILE *wr, uint16_t ival, uint16_t maxval)
  { pnm_check_sample_range(&ival, maxval);
    putc((char)((ival/256)&255), wr);
    putc((char)(ival&255), wr);
  }

uint16_t pnm_quantize(double fval, uint16_t maxval, bool_t isMask, uint32_t badval)
  { /* Maxval with more bits: */
    int N = (uint32_t)maxval;
    int intval; /* Output value with more bits: */
    if (isnan(fval))
      { /* Return {badval} if valid, else {maxval/2}: */
        intval = (badval <= N ? badval : maxval/2);
      }
    else
      { /* Reduce {N} to the effective maxval, with {badval} excluded: */
        if (badval <= N) { N--; }
        /* Convert to integer {intval} in {0..N}: */
        if (fval <= 0.0) 
          { intval = 0; }
        else if (fval >= 1.0)
          { intval = N; }
        else if (isMask)
          { /* Map 0 to 0, 1 to {2*maxval}, linearly, with directed rounding: */ 
            double sval = 2*fval*N;
            intval = (fval < 0.5 ? (int)ceil(sval) : (int)floor(sval) + 1); 
            /* Now round the fraction bit up: */
            intval = intval >> 1;
          }
        else
          { /* Round every interval {[iv/(N+1),(iv+1)/(N+1))} down to {iv}: */
            intval = (int)floor(fval*((double)N + 1));
            /* Should not happen, except perhaps for very large or very small {fdelta}: */
            if (intval > N) { intval = N; }
          }
        assert(intval <= N);
        /* Adjust {intval} for bad-value encoding: */
        if (intval >= badval) { intval++; }
      }
    assert(intval <= maxval);
    return (uint16_t)intval;
  }

double pnm_floatize(uint16_t ival, uint16_t maxval, bool_t isMask, uint32_t badval)
  { pnm_check_sample_range(&ival, maxval);
    uint32_t intval = ival;
    if (intval == badval) { /* Sample is undefined: */ return NAN; }
    /* Effective range, with {badval} excluded: */
    uint32_t N = (uint32_t)maxval;
    if (badval <= N) { N--; }
    assert(N > 0);
    /* Adjust {intval} for bad-value encoding: */
    if (intval > badval) { intval--; }
    assert(intval <= N);
    if (isMask)
      { /* Simple linear scaling: */
        return ((double)intval)/((double)N);
      }
    else
      { /* Map {intval} to the center of its interval: */
        return ((double)intval + 0.5)/((double)N + 1);
      }
  }

double *pnm_make_floatize_table(uint16_t maxval, bool_t isMask, uint32_t badval)
  { uint32_t size = ((uint32_t)maxval) + 1;
    double *cvt = talloc(size, double);
    int intval;
    for (intval = 0; intval <= maxval; intval++)
      { cvt[intval] = pnm_floatize((uint16_t)intval, maxval, isMask, badval); }
    return cvt;
  }

static int pnm_message_count = 0;
  /* Counts calls to {pnm_message}. */

#define PNM_MAX_MESSAGES 100
  /* Abort the program after this many warning messages. */

void pnm_message(char* msg, ...)
  { va_list args;
    va_start(args, msg);
    /* fprintf(stderr, "%s: ", progname); */
    (void)vfprintf(stderr, msg, args);
    fputc('\n', stderr);
    pnm_message_count++;
    if (pnm_message_count > PNM_MAX_MESSAGES) 
      { pnm_error("too many warnings, aborted"); }
    va_end(args);
  }

void pnm_error(char* msg, ...)
  { va_list args;
    va_start(args, msg);
    /* fprintf(stderr, "%s: ", progname); */
    (void)vfprintf(stderr, msg, args);
    fputc('\n', stderr);
    va_end(args);
    exit(1);
  }

bool_t pnm_uint_leq(uint32_t x, uint32_t xmax)
  { return (x <= xmax); }

void pnm_check_sample_range(uint16_t *ival, uint16_t maxval)
  { if (*ival > maxval)
      { pnm_message( "sample value out of bounds (%u > %u)", *ival, maxval);
        *ival = maxval;
      }
  }
            
void pnm_read_pixels
  ( FILE *rd, 
    uint16_t *smp, 
    int cols, 
    int chns,
    uint16_t maxval, 
    bool_t raw,
    bool_t bits
  )
  { assert(pnm_uint_leq(maxval, PNM_FILE_MAX_MAXVAL));
    uint16_t max_byte = 255u;
    uint16_t max_short = 65535u;
    int samples_per_row = cols*chns; /* Samples per row. */
    int k;
    uint16_t *sP;
    if (bits)
      { /* PBM format. */
        /* Note that PBM bit values are inverted -- 0 means white, 1 means black. */
        assert(maxval == 1);
        if (raw)
          { /* Bits packed 8 per byte, big-endianly, complemented, row padded to 8 bits: */
            uint16_t ival;
            for (k = 0, sP = smp; k < samples_per_row; k++)
              { if (k % 8 == 0) { ival = pnm_read_raw_byte(rd, 255); }
                (*sP) = ((ival & 128) == 0); ++sP; 
                ival = (uint16_t)(((int32_t)ival) << 1);
              }
          }
        else
          { /* ASCII digits '0' or '1', with optional whitespace: */
            for (k = 0, sP = smp; k < samples_per_row; k++)
              { uint16_t ival = pnm_read_plain_bit(rd); 
                *sP = ival ^ 1; sP++;
              }
          }
      }
    else 
      { /* PGM/PPM format. */
        if (raw)
          { /* Raw PGM/PPM format (binary samples). */
            if (maxval <= max_byte)
              { /* One byte per sample: */
                for (k = 0, sP = smp; k < samples_per_row; k++)
                  { uint16_t ival = pnm_read_raw_byte(rd, maxval);
                    *sP = ival; ++sP;
                  }
              }
            else if (maxval <= max_short)
              { /* Two bytes per sample: */
                for (k = 0, sP = smp; k < samples_per_row; k++)
                  { uint16_t ival = pnm_read_raw_short(rd, maxval);
                    *sP = ival; ++sP;
                  }
              }
            else
              { pnm_error("maxval (%u) is too large for raw PGM format", maxval); }
          }
        else
          { /* Plain PBM/PGM/PPM format (ASCII decimal samples with whitespace sep). */
            for (k = 0, sP = smp; k < samples_per_row; k++)
              { uint16_t ival = pnm_read_plain_sample(rd, maxval);
                *sP = ival; ++sP;
              }
          }
      }
  }
  
void pnm_write_pixels
  ( FILE* wr, 
    uint16_t *smp, 
    int cols, 
    int chns,
    uint16_t maxval, 
    bool_t raw,
    bool_t bits
  )
  { assert(pnm_uint_leq(maxval, PNM_FILE_MAX_MAXVAL));
    uint16_t max_byte = 255u;
    uint16_t max_short = 65535u;
    int samples_per_row = cols*chns; /* Samples per row. */
    int k;
    uint16_t *sP;
    if (bits)
      { /* PBM format. */
        /* Note that PBM bit values are inverted -- 0 means white, 1 means black. */
        assert(maxval == 1);
        if (raw)
          { /* Raw PBM format - 8 samples per byte, big-endianly, complemented, row padded to byte: */
            uint16_t ival = 0;
            uint16_t msk = 128;
            for (k = 0, sP = smp; k < samples_per_row; k++)
              { assert((*sP) <= 1);
                if ((*sP) == 0) { ival |= msk; }
                msk >>= 1; ++sP;
                if (msk == 0) { pnm_write_raw_byte(wr, ival, 255); msk = 128; ival = 0; }
              }
            /* Flush partially filled byte: */
            if (msk > 0) { pnm_write_raw_byte(wr, ival, 255); }
          }
        else
          { /* Plain PBM format - sample is ASCII '0' or '1', complemented, no whitespace: */
            int maxchars = 70 - 1; /* Max chars in line minus one pixel. */
            int chars = 0; /* Counts characters in current line. */
            for (k = 0, sP = smp; k < samples_per_row; k++)
              { if (chars > maxchars) 
                  { putc('\n', wr); chars = 0; }
                else
                  { putc(' ', wr); chars++; }
                uint16_t ival = *sP; ++sP;
                /* Negate the unit bit only, so that range errors can be detected: */
                pnm_write_plain_bit(wr, ival ^ 1); chars++;
              }
            putc('\n', wr);
          }
      }
    else 
      { /* PGM/PPM format. */
        if (raw)
          { /* Raw PGM/PPM format  (binary pixels). */
            if (maxval <= max_byte)
              { /* Raw single-byte format: */
                for (k = 0, sP = smp; k < samples_per_row; k++)
                  { uint16_t ival = *sP; ++sP;
                    pnm_write_raw_byte(wr, ival, maxval);
                  }
              }
            else if (maxval <= max_short)
              { /* Raw two-byte format: */
                for (k = 0, sP = smp; k < samples_per_row; k++)
                  { uint16_t ival = *sP; ++sP;
                    pnm_write_raw_short(wr, ival, maxval);
                  }
              }
            else
              { pnm_error("maxval (%u) is too large for raw PGM format", maxval); }
          }
        else
          { /* Plain PBM/PGM/PPM format - samples are ASCII decimal with whitespace sep: */
            int maxchars = 70 - 6*chns; /* Max chars in line minus one pixel. */
            int chars = 0; /* Counts characters in current line. */
            for (k = 0, sP = smp; k < samples_per_row; k++)
              { uint16_t ival = *sP; ++sP;
                chars += pnm_write_plain_sample(wr, ival, maxval);
                if (((chars > maxchars) && ((k % chns) == 0)) || (k == samples_per_row-1)) 
                  { putc('\n', wr); chars = 0; }
                else
                  { putc(' ', wr); chars++; }
              }
          }
      }
    fflush(wr);
  }

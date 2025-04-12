/* See fget.h */
/* Last edited on 2025-03-13 06:35:38 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <jsstring.h>
#include <affirm.h>
#include <bool.h>

#include <fget.h>

/* INTERNAL PROTOTYPES */

int32_t fget_digit(FILE *rd, uint32_t base);
  /* If the next character from {rd} exists and is a digit in '0'..'9' or
    a letter in 'a'..'z' or 'A'..'Z' (both meaning 10..35), and its
    numeric value is less than base, consumes that character and returns
    its value. Otherwise leaves the file unchanged and returns -1. */

/* IMPLEMENTATIONS */

bool_t fget_is_formatting_char(char c)
  {
    return 
      (c == ' ') || 
      (c == '\240') || 
      ((c >= '\011') && (c <= '\015'));
  }
  
bool_t fget_is_space(char c)
  {
    return 
      (c == ' ') || 
      (c == '\240') || 
      (c == '\011');
  }

bool_t fget_test_eof(FILE *rd)
  { int32_t r = fgetc(rd);
    if (r == EOF)
      { return TRUE; }
    else
      { ungetc(r, rd);
        return FALSE;
      }
  }

bool_t fget_test_char(FILE *rd, char c)
  { int32_t r = fgetc(rd);
    if (r == EOF)
      { return FALSE; }
    else 
      if (r == (uint8_t)c)
      { return TRUE; }
    else 
      { ungetc(r, rd);
        return FALSE;
      }
  }

void fget_skip_spaces(FILE *rd)
  { int32_t r;
    do 
      { r = fgetc(rd);
        if (r == EOF) { return; }
      }
    while (fget_is_space((char)r));
    ungetc(r, rd); 
    return;
  }

void fget_skip_formatting_chars(FILE *rd)
  { int32_t r;
    do 
      { r = fgetc(rd);
        if (r == EOF) { return; }
      }
    while (fget_is_formatting_char((char)r));
    ungetc(r, rd); 
    return;
  }

void fget_match(FILE *rd, char *t)
  { while ((*t) != '\000')
      { int32_t r = fgetc(rd);
        if ((r == EOF) || ((char)r != (*t)))
          { demand(FALSE, txtcat3("cannot find \"", t, "\"")); }
        t++;
      }
  }

void fget_skip_to_eol(FILE *rd)
  { int32_t r = fgetc(rd);
    while ((r != EOF) && ((char)r != '\n')) { r = fgetc(rd); }
    demand(r != EOF, "no end-of-line at end of file");
    return;
  }

void fget_eol(FILE *rd)
  { int32_t r;
    fget_skip_spaces(rd);
    r = fgetc(rd);
    demand(r != EOF, "no end-of-line at end of file");
    demand((char)r == '\n', "extraneous data on input line");
  }

bool_t fget_test_comment_or_eol(FILE *rd, char cmtc, char **text_P)
  { fget_skip_spaces(rd);
    int32_t r = fgetc(rd);
    if (r == EOF)
      { return FALSE; }
    else if ((char)r == '\n')
      { if (text_P != NULL) { (*text_P) = NULL; }
        return TRUE; 
      }
    else if ((char)r == cmtc)
      { if (text_P != NULL)
          { (*text_P) = fget_line(rd); }
        else 
          { fget_skip_to_eol(rd); }
        return TRUE;
      }
    else
      { ungetc(r, rd);
        return FALSE;
      }
  }

void fget_comment_or_eol(FILE *rd, char cmtc, char **text_P)
  { 
    demand(fget_test_comment_or_eol(rd, cmtc, text_P), "extraneous data on input line");
  }

char fget_char(FILE *rd)
  { fget_skip_spaces(rd);
    int32_t r = fgetc(rd);
    demand(r != EOF, "end of file while looking for item");
    char c = (char)r;
    demand(! fget_is_formatting_char(c), "end of line while looking for item");
    return c;
  }

bool_t fget_bool(FILE *rd)
  { char c = fget_char(rd); 
    if ((c == 't') || (c == 'T') || (c == 'y') || (c == 'Y') || (c == '1')) 
      { return TRUE; }
    else if ((c == 'f') || (c == 'F') || (c == 'n') || (c == 'N') || (c == '0'))
      { return FALSE; }
    else
      { demand(FALSE, "missing or invalid bool_t value"); }
  }

#define INITEXTLENGTH 1024

char *fget_to_delims(FILE *rd, char delim, char *delims)
  { int32_t bufsz = INITEXTLENGTH+1;
    char *buf = talloc(bufsz, char);
    int32_t nb = 0; /* Number of chars stored into {buf}. */
    int32_t r = fgetc(rd);
    while (TRUE)
      { demand(r != EOF, "no end-of-line at end of file");
        char ch = (char)r;
        if (ch == '\n') { break; }
        if (ch == delim) { break; }
        if ((delims != NULL) && (strchr(delims,ch) != NULL)) { break; }
        /* Save char, expanding as needed, and leaving space for final '\0': */
        if (nb >= bufsz-1)
          { bufsz = 2*bufsz;
            buf = retalloc(buf, bufsz, char);
          }
        buf[nb] = ch; nb++;
        r = fgetc(rd);
      }
    assert(r != EOF);
    ungetc(r, rd);
    buf[nb] = '\000'; nb++; 
    if (nb < bufsz) { buf = retalloc(buf, nb, char); }
    return buf;
  }

char *fget_line(FILE *rd)
  { char *s = fget_to_delims(rd, '\n', NULL);
    fget_eol(rd);
    return s;
  }

char *fget_string(FILE *rd)
  { fget_skip_spaces(rd);
    char *s = fget_to_delims(rd, '\n', fget_formatting_chars);
    demand ((*s) != '\000', "item not found");
    return s;
  }

int32_t fget_int32(FILE *rd)
  { int64_t x = fget_int64(rd);
    demand((x >= INT32_MIN) && (x < INT32_MAX), "integer does not fit in {int32_t} type");
    return (int32_t)x;
  }

int64_t fget_int64(FILE *rd)
  { char c = fget_char(rd);
    bool_t positive = TRUE;
    if (c == '+') 
      { /* Ignore. */ }
    else if (c == '-') 
      { positive = FALSE; }
    else
      { ungetc(c, rd); }
    uint64_t x = fget_uint64(rd, 10);
    if (positive) 
      { demand (x <= ((uint64_t)INT64_MAX), "integer does not fit in {int64_t} type");
        return (int64_t)x;
      }
    else
      { demand (x <= ((uint64_t)INT64_MAX) + 1, "integer does not fit in {int64_t} type");
        return (int64_t)(0 - x);
      }
  }

uint32_t fget_uint32(FILE *rd, uint32_t base)
  { uint64_t x = fget_uint64(rd, base);
    demand(x < UINT32_MAX, "integer does not fit in {uint32_t} type");
    return (uint32_t)x;
  }
        
int32_t fget_digit(FILE *rd, uint32_t base)
  { 
    int32_t r = fgetc(rd);
    if (r == EOF) { return -1; }
    int32_t d;
    char c = (char)r;
    if ((c >= '0') && (c <= '9'))
      { d = c - '0'; }
    else if ((c >= 'a') && (c <= 'z'))
      { d = 10 + (c - 'a'); }
    else if ((c >= 'A') && (c <= 'Z'))
      { d = 10 + (c - 'A'); }
    else
      { d = -1;}
    if ((d < 0) || (d >= base))
      { ungetc(r, rd);
        return -1;
      }
    else
      { return d; }
  }

uint64_t fget_uint64(FILE *rd, uint32_t base)
  { demand((base >= 2) && (base <= 36), "invalid base");
    uint64_t xpmax = UINT64_MAX/base; /* Max value of {x} to which a 0 can be appended. */
    fget_skip_spaces(rd);
    int32_t d = fget_digit(rd, base); 
    demand(d >= 0, "invalid number");
    /* Parse digits, append to number: */
    uint64_t x = (uint64_t)d;
    while (TRUE)
      { /* Grab the next digit's value {d}, or {-1} if none: */
        d = fget_digit(rd, base);
        if (d < 0) { break; }
        demand (x <= xpmax, "number does not fit in 64 bits");
        x = base*x;
        demand(x <= UINT64_MAX - (uint64_t)d, "number does not fit in 64 bits");
        x = x + (uint64_t)d;
      }
    return x;
  }

#define MAXNUMLENGTH 1024
  /* Maximum length of a "double" item. */ 

double fget_double(FILE *rd)
  { bool_t debug = FALSE;
    char buf[MAXNUMLENGTH+1];
    int32_t nb = 0; /* Number of chars in {buf}. */
    char c = fget_char(rd);
    /* Parse sign, else put char back: */
    if ((c == '+') || (c == '-'))
      { buf[nb] = c; nb++; 
        if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
      }
    else
      { ungetc(c, rd); }
    /* Now parse number minus the sign. Unset {ok} is missing or bad format. */
    bool_t ok = TRUE; /* Number looks well-formed. */
    int32_t r = fgetc(rd); c = (char)r;
    if (r == EOF)
      { ok = FALSE; }
    else if ((c == 'N') || (c == 'n') || (c == 'I') || (c == 'i'))
      { /* Looks like NaN or INF; gather the whole alpha token: */
        do
          { buf[nb] = c; nb++; 
            if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
            r = fgetc(rd); c = (char)r;
          }
        while 
          ( (r != EOF) && 
            (nb < MAXNUMLENGTH) && 
            (((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z')))
          );
      }
    else if (((c >= '0') && (c <= '9')) || (c == '.'))
      { bool_t has_digs = FALSE;
        /* Assume it is an ordinary number, gobble up int part digits, if any: */
        while ((r != EOF) && (nb < MAXNUMLENGTH) && (c >= '0') && (c <= '9'))
          { buf[nb] = c; nb++; has_digs = TRUE;
            if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
            r = fgetc(rd); c = (char)r; }
        /* Gobble up dot if present: */
        if ((r != EOF) && (nb < MAXNUMLENGTH) && (c == '.'))
          { buf[nb] = c; nb++; 
            r = fgetc(rd); c = (char)r;
            if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
            /* Gobble up fraction part digits, if any: */
            while ((r != EOF) && (nb < MAXNUMLENGTH) && (c >= '0') && (c <= '9'))
              { buf[nb] = c; nb++; has_digs = TRUE; 
                if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
                r = fgetc(rd); c = (char)r;
              }
          }
        if (has_digs)
          { /* Gobble up the exponent code if present: */
            if 
              ( (r != EOF) && 
                has_digs && 
                (nb < MAXNUMLENGTH) && 
                ((c == 'e') || (c == 'E') || (c == 'd') || (c == 'D'))
              )
              { buf[nb] = c; nb++; r = fgetc(rd); c = (char)r;
                /* Gobble up exponent sign if present: */
                if ((r != EOF) && (nb < MAXNUMLENGTH) && ((c == '+') || (c == '-')))
                  { buf[nb] = c; nb++; r = fgetc(rd); c = (char)r; }
                /* Gobble up exponent digits. Assume {strtod} will fail if empty. */
                while ((r != EOF) && (nb < MAXNUMLENGTH) && (c >= '0') && (c <= '9'))
                  { buf[nb] = c; nb++; r = fgetc(rd); c = (char)r; }
              }
          }
        else
          { /* Maybe sign and '.' but no digits: */
            ok = FALSE;
          }
      }
    else
      { /* Possibly a sign, but then a char that can't start a num: */
        ok = FALSE;
      }
    /* Now {r} is {EOF} or the first unused char: */
    if (r != EOF) { ungetc(r, rd); }
    /* Terminate {buf}: */
    buf[nb] = '\000';
    if (debug) { fprintf(stderr, "  buf = \"%s\"\n", buf); }
    /* Parse the number: */
    demand(nb > 0, "expected number, not found"); 
    double x = NAN;
    if (ok)
      { /* Use {strtod} to parse as {double}: */
        char *rest = NULL;
        x = strtod(&(buf[0]), &rest); 
        ok &= ((*rest) == '\000');
      }
    demand(ok, txtcat("invalid number ", &(buf[0])));
    return x;
  }

void fget_skip_spaces_and_match(FILE *rd, char *t)
  { fget_skip_spaces(rd); 
    fget_match(rd, t); 
  }

bool_t fget_skip_spaces_and_test_char(FILE *rd, char c)
  { fget_skip_spaces(rd); 
    return fget_test_char(rd, c); 
  }

void fget_show_next(FILE *wr, char *pref, FILE *rd, char *suff)
  { int32_t chn = fgetc(rd); 
    ungetc(chn, rd);
    fprintf(wr, "%s", pref);
    if (chn == EOF) 
      { fprintf(wr, "EOF"); }
    else if (chn == '\n')
      { fprintf(wr, "LF"); }
    else if (chn == '\r')
      { fprintf(wr, "CR"); }
    else if (chn == '\t')
      { fprintf(wr, "TAB"); }
    else
      { fprintf(wr, "'%c'", chn); }
    fprintf(wr, "%s", suff);
  }

/* Created by J. Stolfi, Unicamp, Dec/2002. */

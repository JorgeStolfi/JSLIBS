/* See fget.h */
/* Last edited on 2023-10-13 13:50:42 by stolfi */

#define _GNU_SOURCE_
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

int32_t fget_digit(FILE *f, uint32_t base);
  /* If the next character from {f} exists and is a digit in '0'..'9' or
    a letter in 'a'..'z' or 'A'..'Z' (both meaning 10..35), and its
    numeric value is less than base, consumes that character and returns
    its value. Otherwise leaves the file unchanged and returns -1. */

/* IMPLEMENTATIONS */

bool_t fget_is_formatting_char(char c)
  {
    return 
      (c == '\000') ||
      (c == ' ') || 
      (c == '\240') || 
      ((c >= '\011') && (c <= '\015'));
  }
  
bool_t fget_is_space(char c)
  {
    return 
      (c == '\000') || 
      (c == ' ') || 
      (c == '\240') || 
      (c == '\011');
  }

bool_t fget_test_eof(FILE *f)
  { int32_t r = fgetc(f);
    if (r == EOF)
      { return TRUE; }
    else
      { ungetc(r, f);
        return FALSE;
      }
  }

bool_t fget_test_char(FILE *f, char c)
  { int32_t r = fgetc(f);
    if (r == EOF)
      { return FALSE; }
    else 
      if (r == (uint8_t)c)
      { return TRUE; }
    else 
      { ungetc(r, f);
        return FALSE;
      }
  }

void fget_skip_spaces(FILE *f)
  { int32_t r;
    do 
      { r = fgetc(f);
        if (r == EOF) { return; }
      }
    while (fget_is_space((char)r));
    ungetc(r, f); 
    return;
  }

void fget_skip_formatting_chars(FILE *f)
  { int32_t r;
    do 
      { r = fgetc(f);
        if (r == EOF) { return; }
      }
    while (fget_is_formatting_char((char)r));
    ungetc(r, f); 
    return;
  }

void fget_match(FILE *f, char *t)
  { while ((*t) != '\000')
      { int32_t r = fgetc(f);
        if ((r == EOF) || ((char)r != (*t)))
          { demand(FALSE, txtcat("cannot find \"", txtcat(t, "\""))); }
        t++;
      }
  }

void fget_skip_to_eol(FILE *f)
  { int32_t r = fgetc(f);
    while ((r != EOF) && ((char)r != '\n')) { r = fgetc(f); }
    demand(r != EOF, "no newline at end of file");
    return;
  }

void fget_eol(FILE *f)
  { int32_t r;
    fget_skip_spaces(f);
    r = fgetc(f);
    demand(r != EOF, "no newline at end of file");
    demand((char)r == '\n', "extraneous data on input line");
  }

void fget_comment_or_eol(FILE *f, char cmtc)
  { fget_skip_spaces(f);
    int32_t r = fgetc(f);
    demand(r != EOF, "unexpected EOF");
    if ((char)r == '\n')
      { return; }
    else if ((char)r == cmtc)
      { do { r = fgetc(f); } while ((r != EOF) && ((char)r != '\n'));
        demand(r != EOF, "no newline at end of file");
      }
    else
      { demand(FALSE, "extraneous data on input line"); }
  }

bool_t fget_test_comment_or_eol(FILE *f, char cmtc)
  { fget_skip_spaces(f);
    int32_t r = fgetc(f);
    if (r == EOF)
      { return FALSE; }
    else if ((char)r == '\n')
      { return TRUE; }
    else if ((char)r == cmtc)
      { do { r = fgetc(f); } while ((r != EOF) && ((char)r != '\n'));
        demand(r != EOF, "no newline at end of file");
        return TRUE;
      }
    else
      { if (r != EOF) { ungetc(r, f); }
        return FALSE;
      } 
  }

char fget_char(FILE *f)
  { fget_skip_spaces(f);
    int32_t r = fgetc(f);
    demand(r != EOF, "end of file while looking for item");
    char c = (char)r;
    demand(! fget_is_formatting_char(c), "end of line while looking for item");
    return c;
  }

bool_t fget_bool(FILE *f)
  { char c = fget_char(f); 
    if ((c == 't') || (c == 'T') || (c == 'y') || (c == 'Y') || (c == '1')) 
      { return TRUE; }
    else if ((c == 'f') || (c == 'F') || (c == 'n') || (c == 'N') || (c == '0'))
      { return FALSE; }
    else
      { demand(FALSE, "missing or invalid bool_t value"); }
  }

#define INITEXTLENGTH 1024

char *fget_to_delim(FILE *f, char del)
  { char dels[2];
    dels[0] = del;
    dels[1] = '\000';
    return fget_to_delims(f, dels);
  }

char *fget_to_delims(FILE *f, char *dels)
  { int32_t bufsz = INITEXTLENGTH+1;
    char *buf = (char *)notnull(malloc(bufsz), "out of mem for buf");;
    int32_t nb = 0; /* Number of chars stored into {buf}. */
    fget_skip_spaces(f);
    int32_t r = fgetc(f);
    while ((r != EOF) && (! fget_is_formatting_char((char)r)) && ((dels == NULL) || (strchr(dels,(char)r) == NULL)))
      { /* Save char, expanding as needed, and leaving space for final '\0': */
        if (nb >= bufsz-1)
          { int32_t bufsz = 2*bufsz;
            buf = (char *)notnull(realloc(buf, bufsz), "out of mem for buf");
          }
        buf[nb] = (char)r; nb++;
        r = fgetc(f);
      }
    if (r != EOF) { ungetc(r, f); }
    buf[nb] = '\000'; nb++; 
    if (nb < bufsz) { buf = (char *)notnull(realloc(buf, nb), "out of mem for result"); }
    return buf;
  }

char *fget_line(FILE *f)
  { char *s = fget_to_delims(f, NULL);
    return s;
  }

char *fget_string(FILE *f)
  { char *s = fget_to_delim(f, '\n');
    demand ((*s) != '\000', "item not found");
    return s;
  }

int32_t fget_int32(FILE *f)
  { int64_t x = fget_int64(f);
    demand((x >= INT32_MIN) && (x < INT32_MAX), "integer does not fit in {int32_t} type");
    return (int32_t)x;
  }

int64_t fget_int64(FILE *f)
  { char c = fget_char(f);
    bool_t positive = TRUE;
    if (c == '+') 
      { /* Ignore. */ }
    else if (c == '-') 
      { positive = FALSE; }
    else
      { ungetc(c, f); }
    uint64_t x = fget_uint64(f, 10);
    if (positive) 
      { demand (x <= ((uint64_t)INT64_MAX), "integer does not fit in {int64_t} type");
        return (int64_t)x;
      }
    else
      { demand (x <= ((uint64_t)INT64_MAX) + 1, "integer does not fit in {int64_t} type");
        return (int64_t)(0 - x);
      }
  }

uint32_t fget_uint32(FILE *f, uint32_t base)
  { uint64_t x = fget_uint64(f, base);
    demand(x < UINT32_MAX, "integer does not fit in {uint32_t} type");
    return (uint32_t)x;
  }
        
int32_t fget_digit(FILE *f, uint32_t base)
  { 
    int32_t r = fgetc(f);
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
      { ungetc(r, f);
        return -1;
      }
    else
      { return d; }
  }

uint64_t fget_uint64(FILE *f, uint32_t base)
  { demand((base >= 2) && (base <= 36), "invalid base");
    uint64_t xpmax = UINT64_MAX/base; /* Max value of {x} to which a 0 can be appended. */
    fget_skip_spaces(f);
    int32_t d = fget_digit(f, base); 
    demand(d >= 0, "invalid number");
    /* Parse digits, append to number: */
    uint64_t x = d;
    while (TRUE)
      { /* Grab the next digit's value {d}, or {-1} if none: */
        d = fget_digit(f, base);
        if (d < 0) { break; }
        demand (x <= xpmax, "number does not fit in 64 bits");
        x = base*x;
        demand(x <= UINT64_MAX - d, "number does not fit in 64 bits");
        x = x + d;
      }
    return x;
  }

#define MAXNUMLENGTH 1024
  /* Maximum length of a "double" item. */ 

double fget_double(FILE *f)
  { bool_t debug = FALSE;
    char buf[MAXNUMLENGTH+1];
    int32_t nb = 0; /* Number of chars in {buf}. */
    char c = fget_char(f);
    /* Parse sign, else put char back: */
    if ((c == '+') || (c == '-'))
      { buf[nb] = c; nb++; 
        if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
      }
    else
      { ungetc(c, f); }
    /* Now parse number minus the sign. Unset {ok} is missing or bad format. */
    bool_t ok = TRUE; /* Number looks well-formed. */
    int32_t r = fgetc(f); c = (char)r;
    if (r == EOF)
      { ok = FALSE; }
    else if ((c == 'N') || (c == 'n') || (c == 'I') || (c == 'i'))
      { /* Looks like NaN or INF; gather the whole alpha token: */
        do
          { buf[nb] = c; nb++; 
            if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
            r = fgetc(f); c = (char)r;
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
            r = fgetc(f); c = (char)r; }
        /* Gobble up dot if present: */
        if ((r != EOF) && (nb < MAXNUMLENGTH) && (c == '.'))
          { buf[nb] = c; nb++; 
            r = fgetc(f); c = (char)r;
            if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
            /* Gobble up fraction part digits, if any: */
            while ((r != EOF) && (nb < MAXNUMLENGTH) && (c >= '0') && (c <= '9'))
              { buf[nb] = c; nb++; has_digs = TRUE; 
                if (debug) { fprintf(stderr, "  parsed '%c'\n", c); }
                r = fgetc(f); c = (char)r;
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
              { buf[nb] = c; nb++; r = fgetc(f); c = (char)r;
                /* Gobble up exponent sign if present: */
                if ((r != EOF) && (nb < MAXNUMLENGTH) && ((c == '+') || (c == '-')))
                  { buf[nb] = c; nb++; r = fgetc(f); c = (char)r; }
                /* Gobble up exponent digits. Assume {strtod} will fail if empty. */
                while ((r != EOF) && (nb < MAXNUMLENGTH) && (c >= '0') && (c <= '9'))
                  { buf[nb] = c; nb++; r = fgetc(f); c = (char)r; }
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
    if (r != EOF) { ungetc(r, f); }
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

void fget_skip_spaces_and_match(FILE *f, char *t)
  { fget_skip_spaces(f); 
    fget_match(f, t); 
  }

bool_t fget_skip_and_test_char(FILE *f, char c)
  { fget_skip_spaces(f); 
    return fget_test_char(f, c); 
  }

/* Created by J. Stolfi, Unicamp, Dec/2002. */

/* See fget.h */
/* Last edited on 2015-11-16 02:35:10 by stolfilocal */

#include <affirm.h>
#include <fget.h>
#include <jsstring.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* INTERNAL PROTOTYPES */

bool_t fget_is_formatting_char(int c);
  /* Returns TRUE iff {c} is a formatting character, namely 
    a space (SPACE, TAB, NUL, NBSP) line break (CR, LF) or page
    break (FF, VT). */
    
bool_t fget_is_space(int c);
  /* Returns TRUE iff {c} is a space char (SPACE, TAB, NUL, NBSP). */

void fget_skip_spaces_to_something(FILE *f);
  /* Skips spaces, then requires that the next character is not 
   a line or page break, or end-of-file. */

int fget_digit(FILE *f, bool_t alpha);
  /* Tries to read the next character from {f}, expecting it to be a
    digit in '0'..'9', or, if {alpha} is true, also a letter in
    'a'..'z' or 'A'..'Z'. If it succeeds, returns the numeric value of
    that character as an integer in 0..9 (if digit) or 10..35 (if
    letter). If the next character is anything else (including EOF),
    puts it back and returns -1. */

/* IMPLEMENTATIONS */

bool_t fget_is_formatting_char(int c)
  {
    return (c == '\000') || (c == ' ') || (c == '\240') || ((c >= '\011') && (c <= '\015'));
  }
  
bool_t fget_is_space(int c)
  {
    return (c == '\000') || (c == ' ') || (c == '\240') || (c == '\011');
  }

bool_t fget_test_char(FILE *f, int c)
  { int r;
    r = fgetc(f);
    if (r == c)
      { return TRUE; }
    else if (r == EOF)
      { return FALSE; }
    else
      { ungetc(r, f); return FALSE; }
  }

void fget_skip_spaces(FILE *f)
  { int c;
    do 
      { c = fgetc(f);
        if (c == EOF) { return; }
      }
    while (fget_is_space(c));
    ungetc(c, f); return;
  }

void fget_skip_formatting_chars(FILE *f)
  { int c;
    do 
      { c = fgetc(f);
        if (c == EOF) { return; }
      }
    while (fget_is_formatting_char(c));
    ungetc(c, f); 
    return;
  }

void fget_skip_to_eol(FILE *f)
  { int c = fgetc(f);
    while ((c != EOF) && (c != '\n')) { c = fgetc(f); }
    demand(c != EOF, "no newline at end of file");
    return;
  }

void fget_match(FILE *f, char *t)
  { while ((*t) != '\000')
      { int c = fgetc(f);
        if (c != (unsigned char)(*t))
          { fprintf(stderr, "next char = '%c'\n", c);
            demand(FALSE, txtcat("cannot find \"", txtcat(t, "\"")));
          }
        t++;
      }
  }

void fget_eol(FILE *f)
  { int c;
    fget_skip_spaces(f);
    c = fgetc(f);
    demand(c != EOF, "no newline at end of file");
    demand(c == '\n', "extraneous data on input line");
  }

void fget_comment_or_eol(FILE *f, int cmtc)
  {
    fget_skip_spaces(f);
    int c = fgetc(f);
    if (c == cmtc)
      { while ((c != '\n') && (c != EOF)) { c = fgetc(f); } }
    demand(c != EOF, "no newline at end of file");
    demand(c == '\n', "extraneous data on input line");
  }

bool_t fget_test_comment_or_eol(FILE *f, int cmtc)
  { 
    fget_skip_spaces(f);
    int c = fgetc(f);
    if (c == cmtc)
      { while ((c != '\n') && (c != EOF)) { c = fgetc(f); } }
    demand(c != EOF, "no newline at end of file");
    if (c == '\n') { return TRUE; }
    ungetc(c, f); 
    return FALSE;
  }

void fget_skip_spaces_to_something(FILE *f)
  { int c;
    fget_skip_spaces(f);
    c = fgetc(f);
    demand((c != EOF) && (! fget_is_formatting_char(c)), "item not found");
    ungetc(c, f);
  }

char fget_char(FILE *f)
  { int c;
    fget_skip_spaces_to_something(f); 
    c = fgetc(f);
    demand(c != EOF, "expecting nonblank char, found end of file");
    return (char)c;
  }

bool_t fget_bool(FILE *f)
  { int c;
    fget_skip_spaces_to_something(f); 
    c = fgetc(f);
    if ((c == 't') || (c == 'T') || (c == 'y') || (c == 'Y') || (c == '1')) 
      { return TRUE; }
    else if ((c == 'f') || (c == 'F') || (c == 'n') || (c == 'N') || (c == '0'))
      { return FALSE; }
    else
      { ungetc(c, f);
        demand(FALSE, "missing or invalid bool_t value");
        return FALSE;
      }
  }

#define INITEXTLENGTH 1024

char *fget_to_delim(FILE *f, char *dels)
  { char buf[INITEXTLENGTH+1];
    char *pbuf = &(buf[0]);
    int bufsz = INITEXTLENGTH+1;
    int i = 0, c;
    fget_skip_spaces(f);
    c = fgetc(f);
    while ((c != EOF) && (! fget_is_formatting_char(c)) && ((dels == NULL) || (strchr(dels,c) == NULL)))
      { if (i >= bufsz-1)
          { int tmpsz = 2*bufsz;
            char *tmp = (char *)notnull(malloc(tmpsz), "out of mem for buf");
            pbuf[i] = '\000';
            strcpy(tmp, pbuf);
            if (pbuf != &(buf[0])) { free(pbuf); }
            pbuf = tmp; bufsz = tmpsz;
          }
        pbuf[i] = (char)c; i++;
        c = fgetc(f);
      }
    if (c != EOF) { ungetc(c, f); }
    pbuf[i] = '\000'; i++; 
    if ((i < bufsz) || (pbuf == &(buf[0])))
      { char *tmp = (char *)notnull(malloc(i), "out of mem for result");
        strcpy(tmp, pbuf);
        if (pbuf != &(buf[0])) { free(pbuf); }
        pbuf = tmp;
      }
    return pbuf;
  }

char *fget_string(FILE *f)
  { char *s = fget_to_delim(f, NULL);
    demand ((*s) != '\000', "item not found");
    return s;
  }

int fget_digit(FILE *f, bool_t alpha)
  { int c = fgetc(f);
    if ((c >= '0') && (c <= '9'))
      { return c - '0'; }
    else if (alpha)
      { if ((c >= 'a') && (c <= 'z'))
          { return 10 + (c - 'a'); }
        else if ((c >= 'A') && (c <= 'Z'))
          { return 10 + (c - 'A'); }
      }
    /* No digit found: */
    if (c != EOF) { ungetc(c, f); } 
    return -1;
  }

int fget_int(FILE *f)
  { int64_t x = fget_int64(f);
    demand((x >= INT_MIN) && (x < INT_MAX), "integer does not fit in {int} type");
    return (int)x;
  }

int fget_int32(FILE *f)
  { int64_t x = fget_int64(f);
    demand((x >= INT32_MIN) && (x < INT32_MAX), "integer does not fit in {int32_t} type");
    return (int)x;
  }

int64_t fget_int64(FILE *f)
  { int c;
    bool_t positive = TRUE;
    fget_skip_spaces_to_something(f);
    c = fgetc(f);
    if (c == '+') 
      { c = fgetc(f); }
    else if (c == '-') 
      { positive = FALSE; c = fgetc(f); }
    if (c != EOF) { ungetc(c, f); }
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

unsigned int fget_uint(FILE *f, int base)
  { uint64_t x = fget_uint64(f, base);
    demand(x < UINT_MAX, "integer does not fit in {unsigned int} type");
    return (unsigned int)x;
  }

unsigned int fget_uint32(FILE *f, int base)
  { uint64_t x = fget_uint64(f, base);
    demand(x < UINT32_MAX, "integer does not fit in {uint32_t} type");
    return (unsigned int)x;
  }

uint64_t fget_uint64(FILE *f, int base)
  { demand((base >= 2) && (base <= 36), "invalid base");
    uint64_t x = 0;
    fget_skip_spaces_to_something(f);
    int d = fget_digit(f, (base > 10));
    demand(d >= 0, "number not found"); 
    do
      { demand(d < base, "invalid digit in number"); 
        demand (x <= (UINT64_MAX - d)/base, "number does not fit in 64 bits");
        x = base*x + d;
        d = fget_digit(f, (base > 10));
      }
    while (d >= 0);
    return x;
  }

#define MAXNUMLENGTH 1024
  /* Maximum length of a "double" item. */ 

double fget_double(FILE *f)
  { char buf[MAXNUMLENGTH+1];
    int i, c;
    double x;
    char *rest;
    fget_skip_spaces_to_something(f);
    i = 0;
    c = fgetc(f);
    if ((c == '+') || (c == '-'))
      { buf[i] = (char)c; i++; c = fgetc(f); }
    if ((c != EOF) && ((c == 'N') || (c == 'n') || (c == 'I') || (c == 'i')))
      { /* Looks like NaN or INF; gather the whole alpha token: */
        buf[i] = (char)c; i++; c = fgetc(f); 
        while ((i < MAXNUMLENGTH) && (c != EOF) && (((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z'))))
          { buf[i] = (char)c; i++; c = fgetc(f); }
      }
    else
      { /* Assume it is an ordinary number: */
        while ((i < MAXNUMLENGTH) && (c != EOF) && (c >= '0') && (c <= '9'))
          { buf[i] = (char)c; i++; c = fgetc(f); }
        if ((i < MAXNUMLENGTH) && (c == '.'))
          { buf[i] = (char)c; i++; c = fgetc(f); }
        while ((i < MAXNUMLENGTH) && (c >= '0') && (c <= '9'))
          { buf[i] = (char)c; i++; c = fgetc(f); }
        if ((i < MAXNUMLENGTH) && ((c == 'e') || (c == 'E') || (c == 'd') || (c == 'D')))
          { buf[i] = (char)c; i++; c = fgetc(f);
            if ((i < MAXNUMLENGTH) && ((c == '+') || (c == '-')))
              { buf[i] = (char)c; i++; c = fgetc(f); }
            while ((i < MAXNUMLENGTH) && (c >= '0') && (c <= '9'))
              { buf[i] = (char)c; i++; c = fgetc(f); }
          }
      }
    /* Now {c} is {EOF} or the first unused char: */
    if (c != EOF) { ungetc(c, f); }
    buf[i] = '\000'; 
    x = strtod(&(buf[0]), &rest);
    demand((*rest) == '\000', txtcat("invalid number", &(buf[0])));
    return x;
  }

void fget_skip_and_match(FILE *f, char *t)
  { fget_skip_spaces(f); 
    fget_match(f, t); 
  }

bool_t fget_skip_and_test_char(FILE *f, int c)
  { fget_skip_spaces(f); 
    return fget_test_char(f, c); 
  }

/* Created by J. Stolfi, Unicamp, Dec/2002. */

/* See jsstring.h */
/* Last edited on 2023-11-25 10:30:01 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <jsstring.h>

int32_t isprefix(const char *s, const char *t)
  { while (((*s)!='\000') &&((*t)!='\000') && ((*s) == (*t)))
      { s++; t++; }
    return (*s) == '\000';
  }

char *prefix(const char *s, int32_t len)
  { demand((len >= 0) && (len <= strlen(s)), "invalid {len}");
    char *r = NULL;
    asprintf(&r, "%.*s", len, s);
    assert((r != NULL) && strlen(r) == len);
    return r;
  }

char *txtcat (const char *a, const char *b)
  { char *r = NULL;
    asprintf(&r, "%s%s", a, b);
    return (char*)notnull(r, "no mem");
  }

char *txtcat3 (const char *a, const char *b, const char *c)
  { char *r = NULL;
    asprintf(&r, "%s%s%s", a, b, c);
    return (char*)notnull(r, "no mem");
  }

char *txtcat4 (const char *a, const char *b, const char *c, const char *d)
  { char *r = NULL;
    asprintf(&r, "%s%s%s%s", a, b, c, d);
    return (char*)notnull(r, "no mem");
  }
  
char *txtrep(const char* x, uint32_t n)
  { uint64_t m = strlen(x);
    char *r = talloc(n*m  +  1, char);
    char *p = r;
    uint32_t k;
    for (k = 0; k < n; k++) { strcpy(p, x); p += m; }
    return r;
  }  

char *add_ext(const char *name, const char *ext)
  { 
    if ((strcmp(name, "") == 0) || (strcmp(name, "-") == 0))
      { return txtcat(name, ""); }
    else
      { return txtcat(name, ext); }
  }
  
#define is_space(ch) (((ch) == ' ') || ((ch) == '\240') || ((ch) == '\011'))

char *trim_spaces(char *x, bool_t at_beg, bool_t at_end)
  { uint64_t m = strlen(x);
    char *p = x;
    if (at_beg) { while (is_space(*p)) { p++; } }
    char *q = x + m;
    if (at_end) { while ((q > p) && is_space(*(q-1))) { q--; } }
    assert(p <= q);
    uint64_t n = (uint64_t)(q - p);
    char *r = talloc(n+1, char);
    bcopy(p, r, n); 
    *(r + n) = '\000';
    return r;
  }

char *fmt_int(int64_t x, uint32_t wid)
  { char *r = NULL;
    asprintf(&r, "%0*ld", wid, x);
    return r;
  }

char *escapify(char *x)
  { uint64_t m = strlen(x);
    char *r = talloc(4*m + 1, char);
    char *px = x;
    char *pr = r;
    while ((*px) != 0)
      { char ch = (*px); 
        if (ch == '\\')
          { (*pr) = '\\'; pr++;
            (*pr) = '\\'; pr++;
          }
        else if (ch == '\012')
          { (*pr) = '\\'; pr++;
            (*pr) = 'n'; pr++;
          }
        else if (ch == '\011')
          { (*pr) = '\\'; pr++;
            (*pr) = 't'; pr++;
          }
        else if ((ch >= ' ') && (ch <= '~'))
          { (*pr) = ch; pr++; }
        else
          { (*pr) = '\\'; pr++;
            (*pr) = (char)('0' + (((uint8_t)ch) / 64)); pr++;
            (*pr) = (char)('0' + ((((uint8_t)ch) / 8) % 8));; pr++;
            (*pr) = (char)('0' + (((uint8_t)ch) % 8));; pr++;
          }
        px++;
      }
    (*pr) = 0; pr++;
    r = realloc(r, pr - r);
    return r;
  }


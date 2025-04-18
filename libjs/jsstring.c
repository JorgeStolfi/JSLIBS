/* See jsstring.h */
/* Last edited on 2025-03-13 08:58:19 by stolfi */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <jsprintf.h>
#include <bool.h>

#include <jsstring.h>

int32_t isprefix(const string_t s, const string_t t)
  { char *p = s;
    char *q = t;
    while (((*p)!='\000') &&((*q)!='\000') && ((*p) == (*q)))
      { p++; q++; }
    return (*p) == '\000';
  }

string_t prefix(const string_t s, int32_t len)
  { demand((len >= 0) && (len <= strlen(s)), "invalid {len}");
    string_t r = jsprintf("%.*s", len, s);
    assert((r != NULL) && strlen(r) == len);
    return r;
  }

string_t txtcat (const string_t a, const string_t b)
  { char *aa = (a == NULL ? "" : a);
    char *bb = (b == NULL ? "" : b);
    string_t r = jsprintf("%s%s", aa, bb);
    return (string_t)notnull(r, "no mem");
  }

string_t txtcat3 (const string_t a, const string_t b, const string_t c)
  { char *aa = (a == NULL ? "" : a);
    char *bb = (b == NULL ? "" : b);
    char *cc = (c == NULL ? "" : c);
    string_t r = jsprintf("%s%s%s", aa, bb, cc);
    return (string_t)notnull(r, "no mem");
  }

string_t txtcat4 (const string_t a, const string_t b, const string_t c, const string_t d)
  { char *aa = (a == NULL ? "" : a);
    char *bb = (b == NULL ? "" : b);
    char *cc = (c == NULL ? "" : c);
    char *dd = (d == NULL ? "" : d);
    string_t r = jsprintf("%s%s%s%s", aa, bb, cc, dd);
    return (string_t)notnull(r, "no mem");
  }
  
string_t txtrep(const string_t  x, uint32_t n)
  { char *xx = (x == NULL ? "" : x);
    uint64_t m = strlen(xx);
    string_t r = talloc(n*m  +  1, char);
    string_t p = r;
    uint32_t k;
    for (k = 0; k < n; k++) { strcpy(p, xx); p += m; }
    return r;
  }  

string_t add_ext(const string_t name, const string_t ext)
  { 
    if ((strcmp(name, "") == 0) || (strcmp(name, "-") == 0))
      { return txtcat(name, ""); }
    else
      { return txtcat(name, ext); }
  }
  
#define is_space(ch) (((ch) == ' ') || ((ch) == '\240') || ((ch) == '\011'))

string_t trim_spaces(string_t x, bool_t at_beg, bool_t at_end)
  { uint64_t m = strlen(x);
    char *p = x;
    if (at_beg) { while (is_space(*p)) { p++; } }
    char *q = x + m;
    if (at_end) { while ((q > p) && is_space(*(q-1))) { q--; } }
    assert(p <= q);
    uint64_t n = (uint64_t)(q - p);
    string_t r = talloc(n+1, char);
    bcopy(p, r, n); 
    *(r + n) = '\000';
    return r;
  }

string_t fmt_int(int64_t x, uint32_t wid)
  { string_t r = jsprintf("%0*ld", wid, x);
    return r;
  }

string_t escapify(string_t x)
  { uint64_t m = strlen(x);
    string_t r = talloc(4*m + 1, char);
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
    r = realloc(r, (size_t)(pr - r));
    return r;
  }


/* See jsstring.h */
/* Last edited on 2019-12-14 17:01:59 by jstolfi */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <jsstring.h>

#include <affirm.h>

int32_t isprefix(const char *s, const char *t)
  { while (((*s)!='\000') &&((*t)!='\000') && ((*s) == (*t)))
      { s++; t++; }
    return (*s) == '\000';
  }

char *txtcat (const char *a, const char *b)
  { char *r = malloc(strlen(a) + strlen(b) + 1);
    affirm (r != NULL, "memory exhausted");
    strcpy(r, a);
    strcat(r, b);
    return(r);
  }

char *txtcat3 (const char *a, const char *b, const char *c)
  { char *r = malloc(strlen(a) + strlen(b) + strlen(c) + 1);
    affirm (r != NULL, "memory exhausted");
    strcpy(r, a);
    strcat(r, b);
    strcat(r, c);
    return(r);
  }
  
char *txtrep(const char* x, uint32_t n)
  { uint64_t m = strlen(x);
    char *r = malloc(n*m  +  1);
    char *p = r;
    uint32_t k;
    for (k = 0; k < n; k++) { strcpy(p, x); p += m; }
    return r;
  }  

char *addext(const char *name, const char *ext)
  { 
    if ((strcmp(name, "") == 0) || (strcmp(name, "-") == 0))
      { return txtcat(name, ""); }
    else
      { return txtcat(name, ext); }
  }

#define fmtint_BUFSIZE 1000
char *fmtint(int64_t x, uint32_t wid)
  {
    demand(wid <= fmtint_BUFSIZE-1, "{wid} too big");
    char buf[fmtint_BUFSIZE];
    int32_t rcode = snprintf(buf, fmtint_BUFSIZE, "%0*ld", wid, x);
    assert(rcode >= 0);
    { uint64_t n = strlen(buf);
      char *res = (char *)malloc(n + 1);
      strcpy(res, buf);
      return res;
    }
#undef fmtint_BUFSIZE   
  }



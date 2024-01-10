/* See jsstring.h */
/* Last edited on 2013-10-25 01:08:09 by stolfilocal */

#include <jsstring.h>

#include <affirm.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int isprefix(const char *s, const char *t)
  { while (((*s)!='\000') &&((*t)!='\000') && ((*s) == (*t)))
      { s++; t++; }
    return (*s) == '\000';
  }

char *txtcat (const char *a, const char *b)
  { char *r = malloc(strlen(a)+strlen(b)+1);
    affirm (r != NULL, "memory exhausted");
    strcpy(r, a);
    strcat(r, b);
    return(r);
  }

char *txtcat3 (const char *a, const char *b, const char *c)
  { char *r = malloc(strlen(a)+strlen(b)+strlen(c)+1);
    affirm (r != NULL, "memory exhausted");
    strcpy(r, a);
    strcat(r, b);
    strcat(r, c);
    return(r);
  }

char *addext(const char *name, const char *ext)
  { 
    if ((strcmp(name, "") == 0) || (strcmp(name, "-") == 0))
      { return txtcat(name, ""); }
    else
      { return txtcat(name, ext); }
  }

char *fmtint(int x, unsigned int wid)
  {
#define fmtint_BUFSIZE 200
    char buf[fmtint_BUFSIZE];
    int rcode = snprintf(buf, fmtint_BUFSIZE, "%0*d", wid, x);
    affirm (rcode >= 0, "snprintf failed");
    { int n = (int)strlen(buf);
      char *res = (char *)malloc(n+1);
      strcpy(res, buf);
      return res;
    }
#undef fmtint_BUFSIZE   
  }



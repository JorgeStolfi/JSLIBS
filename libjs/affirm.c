/* See affirm.h */
/* Last edited on 2024-11-20 08:19:12 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <affirm.h>

void programerror (const char *msg, const char *file, int32_t line, const char* proc)
  { fprintf (stderr, "%s:%d: ** (%s) %s\n", file, line, proc, msg);
    exit(1);
  }

void *checknotnull(void *p, const char *msg, const char *file, int32_t line, const char* proc)
  { if (p == NULL) { programerror(msg, file, line, proc); }
    return (void *)p;
  }

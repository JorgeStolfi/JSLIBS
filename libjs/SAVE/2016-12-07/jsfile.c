/* See jsfile.h */
/* Last edited on 2014-06-09 17:31:31 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <jsfile.h>
#include <affirm.h>

FILE *open_write(const char *name, bool_t verbose)
  {
    FILE *f;
    demand (name != NULL, "output file name is null");
    demand (strlen(name) != 0, "output file name is empty");
    if (strcmp(name, "-") == 0)
      { if (verbose) { fprintf(stderr, "writing to standard output\n"); fflush(stderr); }
        f = stdout; 
      }
    else
      { if (verbose) { fprintf(stderr, "writing to %s\n", name); fflush(stderr); }
        f = fopen(name, "w");
      }
    affirm (f != NULL, "could not open output file");
    return f;
  }

FILE *open_read(const char *name, bool_t verbose)
  {
    FILE *f;
    demand (name != NULL, "input file name is null");
    demand (strlen(name) != 0, "input file name is empty");
    if (strcmp(name, "-") == 0)
      { if (verbose) { fprintf(stderr, "reading from standard input\n"); fflush(stderr); }
        f = stdin;
      }
    else
      { if (verbose) { fprintf(stderr, "reading from %s\n", name); fflush(stderr); }
        f = fopen(name, "r");
      }
    affirm (f != NULL, "could not open input file");
    return f;
  }

FILE *open_read_tag_ext(char *name, char *tag, char *ext, bool_t verbose)
  { 
    char *fileName = NULL;
    asprintf(&fileName, "%s%s%s", name, tag, ext);
    FILE *rd = open_read(fileName, verbose);
    free(fileName);
    return rd;
  }

FILE *open_write_tag_ext(char *name, char *tag, char *ext, bool_t verbose)
  { 
    char *fileName = NULL;
    asprintf(&fileName, "%s%s%s", name, tag, ext);
    FILE *wr = open_write(fileName, verbose);
    free(fileName);
    return wr;
  }

char *read_line(FILE *f)
  {
    int mc = 0;
    int nc = 0;
    char *s = NULL;
    int c;
    do
      { c = getc(f); 
        if ((c == EOF) && (nc == 0)) return(NULL);
        if ((c == EOF) || (c == '\n')) c = '\000';
        if (c == '\t') c = ' ';
	if (nc >= mc)
          { if (mc == 0)
	      { mc = 40; s = (char *) malloc(mc*sizeof(char)); }
	    else 
	      { mc *= 2; s = (char *) realloc ((void *) s, mc*sizeof(char)); }
          }
	affirm (s != NULL, "alloc failed");
	s[nc] = (char)c;
	nc++;
      }
    while (c != '\000');
    return (s);
  }


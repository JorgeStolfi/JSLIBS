/* See filefmt.h */
/* Last edited on 2015-11-16 02:35:37 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>

#include <filefmt.h>
#include <fget.h>
#include <jsstring.h>
#include <affirm.h>
#include <vec.h>

void filefmt_write_header(FILE *wr, char *type, char *version)
  { fprintf(wr, "begin %s",  type);
    if (version != NULL) { fprintf(wr, " (format of %s)",  version); }
    fprintf(wr, "\n");
  }

void filefmt_write_footer(FILE *wr, char *type)
  { fprintf(wr, "end %s\n",  type);
  }

char *filefmt_make_header(char *type, char *version) 
  { char *h = NULL;
    if (version != NULL)
      { asprintf(&h, "begin %s (format of %s)\n", type, version); }
    else
      { asprintf(&h, "begin %s\n", type); }
    return h;
  }

char *filefmt_make_footer(char *type)
  { char *h = NULL;
    asprintf(&h, "end %s\n", type);
    return h;
  }

void filefmt_read_gen_header(FILE *rd, char **typeP, char **versionP)
  {
    /* Parse "begin {TYPE}": */
    fget_skip_formatting_chars(rd);
    fget_match(rd, "begin");
    fget_skip_spaces(rd);
    char *type = fget_string(rd);
    fget_skip_spaces(rd);
    
    /* Parse optional "(format of {VERSION})": */
    char *version = NULL;
    if (fget_test_char(rd, '('))
      { fget_match(rd, "format of");
        fget_skip_spaces(rd);
        version = fget_to_delim(rd, ")");
        fget_skip_spaces(rd);
        fget_match(rd, ")");
        fget_skip_spaces(rd);
      }
    
    /* Parses end-of-line: */
    fget_match(rd, "\n");
    
    /* Return what we found: */
    (*typeP) = type;
    (*versionP) = version;
  }
  
void filefmt_read_gen_footer(FILE *rd, char **typeP)
  { 
    fget_skip_formatting_chars(rd);
    fget_match(rd, "end ");
    fget_skip_spaces(rd);
    char *type = fget_string(rd);
    fget_skip_spaces(rd);
    fget_match(rd, "\n");
    
    (*typeP) = type;
  }

void filefmt_read_header(FILE *rd, char *type, char *version)
  { 
    char *rd_type = NULL;
    char *rd_version = NULL;
    filefmt_read_gen_header(rd, &rd_type, &rd_version);
    demand(strcmp(type, rd_type) == 0, "wrong file type in header");
    if (version != NULL)
      { demand(strcmp(version, rd_version) == 0, "wrong file version in header"); }
    else
      { demand(rd_version == NULL, "spurious file version in header"); }
    free(rd_type);
    free(rd_version);
  }

void filefmt_read_footer(FILE *rd, char *type)
  { 
    char *rd_type = NULL;
    filefmt_read_gen_footer(rd, &rd_type);
    demand(strcmp(type, rd_type) == 0, "wrong file type in footer");
    free(rd_type);
  }

void filefmt_write_comment(FILE *wr, char *cmt, char prefix)
  { if (cmt == NULL) { return; }
    while(*cmt != '\000')
      { /* Write a new line, advance {cmt} to start of next one or to '\000' */
        char c = (*cmt);
        fputc(prefix, wr);
        fputc(' ', wr);
        fputc(c, wr);
        while (c != '\n')
          { cmt++; c = (*cmt); 
            if (c == '\000') { c = '\n'; }
            fputc(c, wr);
          }
        if (*cmt == '\n') { cmt++; }
      }
  }

char *filefmt_read_comment(FILE *rd, char prefix)
  { char_vec_t cmt = char_vec_new(200);
    int ncmt = 0; /* Number of characters read */
    while (TRUE)
      { int c;
        /* Try to get the prefix char: */
        c = fgetc(rd); 
        if (c == EOF) { break; }
        if (c != prefix) { ungetc(c, rd); break; }
        /* Skip optional blank after prefix char: */
        c = fgetc(rd); 
        demand(c != EOF, "missing final newline in comment");
        if (c != ' ') { ungetc(c, rd); }
        /* Read comment text until end-of-line: */
        do
          { c = fgetc(rd); 
            demand(c != EOF, "missing final newline in comment");
            char_vec_expand(&cmt, ncmt);
            cmt.e[ncmt] = (char)c;
            ncmt ++;
          }
        while (c != '\n');
      }
    char_vec_expand(&cmt, ncmt);
    cmt.e[ncmt] = '\000';
    ncmt ++;
    char_vec_trim(&cmt, ncmt);
    return cmt.e;
  }

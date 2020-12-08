/* See {stpoly_STP.h} */
/* Last edited on 2016-04-18 19:53:34 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bool.h>
#include <argparser.h>
#include <fget.h>
#include <jsmath.h>

#include <stpoly_STP.h>

bool_t stpoly_STP_is_break(int ch);
  /* TRUE iff {ch} is a line or page break character: ASCII CR ('\015'), 
    LF ('\012'), VT ('\013'), or FF ('\014'). */

bool_t stpoly_STP_is_space(int ch);
  /* TRUE iff {ch} is a blank space character, but not a line or page
    break; namely, if {ch} is ASCII SP ('\040'), TAB ('\011'), NUL
    ('\000') and ISO Latin non-breaking space ('\240'). */
    
/* STP FILE PARSING */

/*
  The following procedures  increment {*lineP} whenever they skip over a line or page break character
  (except that Apple's CR-LF pair is counted as a single line break).  The {fileName}
  argument is used in error messages. */

bool_t stpoly_STP_skip_white_space(FILE *rd, char *fileName, int *lineP);
  /* Skips space characters, line breaks, and page breaks from {rd}.  Returns TRUE if it finds
    a character other than those (which is not consumed), FALSE if it runs into
    the end-of-file.  */

char *stpoly_STP_get_keyword(FILE *rd, char *fileName, int *lineP);
  /* Skips space characters, line breaks, and page breaks from {rd}. If
    it runs into end-of-file or a character that is not an ASCII letter,
    fails with an error message. Otherwise reads characters from {rd},
    until end-of-file, space, line break, or page break, and returns
    those characters as a newly allocated string. */

void stpoly_STP_check_keyword(FILE *rd, char *fileName, int *lineP, char *key);
  /* Skips space characters, line breaks, and page breaks from {rd}
    and obtains the next token with {stpoly_STP_get_keyword}.
    Checks whether it is equal to {key}. Fails with an error message if
    there is no next token in {rd}, or the next token does not start
    with an ASCII letter, or it is not equal to {key}. */

bool_t stpoly_STP_ascii_read_edge(FILE *rd, char *fileName, int *lineP, stpoly_STP_edge_t *edge);
  /* Tries to read another line segment from the ASCII STP file {rd}. If the
    next token in {rd} is "side", reads the segment data up to and
    including the "endside" token, stores that data in {*edge}, and
    returns {TRUE}. If the next token is "endpoly" instead, consumes it
    and returns {FALSE}. Skips spaces, line breaks, and page breaks, and
    other formatting chars before and inside the parsed data. Fails with
    error message if the next token is something else, or the edge
    is malformed. */

int stpoly_STP_binary_read_header(FILE *rd, char *fileName, int *lineP);
  /* Reads the 80-byte header of the binary STP file {rd}, and the 
    next line with the number of edges {ne}.  Returns {ne}.
    The line counter {*lineP} is incremented by 2. */

void stpoly_STP_binary_read_edge(FILE *rd, char *fileName, int *lineP, stpoly_STP_edge_t *edge);
  /* Tries to read another line segment from the binary STP file {rd}. 
    Assumes that each edge is 16 bytes: two endpoints, each with a {X} and a {Y}
    coordinate (4*float). Increments the line counter {*lineP} by 1. */

/* IMPLEMENTATIONS */

bool_t stpoly_STP_is_break(int ch)
  { 
    return (ch == '\012') || (ch == '\013') || (ch == '\014') || (ch == '\015');
  }

bool_t stpoly_STP_is_space(int ch)
  { 
    return (ch == '\000') || (ch == '\011') || (ch == ' ') || (ch == '\240');
  }

bool_t stpoly_STP_skip_white_space(FILE *rd, char *fileName, int *lineP)
  {
    int prev_ch = EOF;  /* Previous character parsed in this call, or EOF if first. */
    int ch = fgetc(rd);
    while (ch != EOF)
      { if (stpoly_STP_is_break(ch))
          { /* Possible line break, unless it is part of an Apple end-of-line (LF/VT/FF after a CR): */
            if ((prev_ch != '\015') || (ch != '\012')) { (*lineP)++; }
          }
        else if (! stpoly_STP_is_space(ch))
          { ungetc(ch, rd); return TRUE; }
        ch = fgetc(rd);
      }
    return FALSE;
  }
  
char *stpoly_STP_get_keyword(FILE *rd, char *fileName, int *lineP)
  {
    if (! stpoly_STP_skip_white_space(rd, fileName, lineP))
      { fprintf(stderr, "%s:%d: ** expecting keyword, found end-of-file\n", fileName, (*lineP)); 
        exit(1);
      }
    int ch = fgetc(rd);
    if (((ch < 'a') || (ch > 'z')) && ((ch < 'A') || (ch > 'Z')))
      { fprintf(stderr, "%s:%d: ** expecting keyword, found '%c'\n", fileName, (*lineP), ch); 
        exit(1);
      }
    ungetc(ch, rd);
    return fget_string(rd);
  }

void stpoly_STP_check_keyword(FILE *rd, char *fileName, int *lineP, char *key)
  {
    char *tok = stpoly_STP_get_keyword(rd, fileName, lineP);
    if (strcmp(tok, key) != 0)
      { fprintf(stderr, "%s:%d: ** expecting '%s', found '%s'\n", fileName, (*lineP), key, tok);
        exit(1);
      }
    free(tok);
  }

void stpoly_STP_read(char *fileName, bool_t binary, stpoly_STP_edge_proc_t *process_edge)  
  {
    char *ftype = (binary ? "rb" : "r");
    FILE *rd = fopen(fileName, ftype);
    if (rd == NULL) 
      { fprintf(stderr, "** failed to open file '%s'\n", fileName);
        exit(1);
      }
    stpoly_STP_edge_t edge;
    int line = 1; /* Line number in file. */
    int ne = 0;   /* Number of segs read from the STP file. */
    
    if (binary)
      { int ne = stpoly_STP_binary_read_header(rd, fileName, &line);
        int it;
        for (it = 0; it < ne; it++) 
          { stpoly_STP_binary_read_edge(rd, fileName, &line, &edge);
            process_edge(line, &edge);
          }
      }
    else
      { stpoly_STP_check_keyword(rd, fileName, &line, "poly");
        while (stpoly_STP_ascii_read_edge(rd,fileName, &line, &edge))
          { ne++;
            process_edge(line, &edge);
          }
      }
  }

bool_t stpoly_STP_ascii_read_edge(FILE *rd, char *fileName, int *lineP, stpoly_STP_edge_t *edge)
  { 
    bool_t debug = TRUE;
    
    char *tok = stpoly_STP_get_keyword(rd, fileName, lineP);
    if (debug) { fprintf(stderr, "tok = %s\n", tok); }
    
    if (strcmp(tok, "endpoly") == 0)
      { return FALSE; }
    else if (strcmp(tok, "side") == 0)
      { int i, k;
        for (k = 0; k < 2; k++) 
          { stpoly_STP_check_keyword(rd, fileName, lineP, "vertex"); 
            for (i = 0; i < 2; i++) 
              { if (! stpoly_STP_skip_white_space(rd, fileName, lineP))
                  { fprintf(stderr, "%s:%d: ** expecting vertex coord, found end-of-file\n", fileName, (*lineP));
                    exit(1);
                  }
                edge->v[k].c[i] = (float)fget_double(rd);
              }
          }
        stpoly_STP_check_keyword(rd, fileName, lineP, "endside");
        return TRUE;
      }
    else
      { fprintf(stderr, "%s:%d: ** expected 'side' or 'endpoly', found '%s'\n", fileName, (*lineP), tok);
        exit(1);
      }
    fclose(rd);
  }
  
int stpoly_STP_binary_read_header(FILE *rd, char *fileName, int *lineP)
  { 
    char title[80];
    int ne;
    int err;
    err = (int)fread(title, 80, 1, rd);
    if (err != 1) { fprintf(stderr, "%s: error reading binary STP file header\n", fileName); exit(1); }
    (*lineP)++;
    err = (int)fread((void*)(&ne), 4, 1, rd);
    if (err != 1) { fprintf(stderr, "%s: error reading num edges from binary STP file\n", fileName); exit(1); }
    if ((ne < 0) || (ne > stpoly_ne_MAX))
      { fprintf(stderr, "%s:%d: ** invalid number of edges (%d) in binary file\n", fileName, (*lineP), ne);
        exit(1);
      }
    (*lineP)++;
    return ne;
  }

void stpoly_STP_binary_read_edge(FILE *rd, char *fileName, int *lineP, stpoly_STP_edge_t *edge)
  { 
    /* Read the coordinates of 2 vertices (2*2 floats) into {vc[0..3]}: */
    float vc[4];
    int err;
    int k;
    for (k = 0; k < 4; k++) 
      { err = (int)fread((void*)(&vc[k]), sizeof(float), 1, rd);
        if (err != 1) { fprintf(stderr, "%s:%d: error reading binary STP datum %d\n", fileName, (*lineP), k); exit(1); }
      }
    
    /* Fill the {stpoly_STP_edge_t}. */
    edge->v[0] = (r2_t){{  (double)vc[0],  (double)vc[1] }};
    edge->v[1] = (r2_t){{  (double)vc[2],  (double)vc[3] }};
  }

void stpoly_STP_print_edge(FILE *wr, stpoly_STP_edge_t *edge)
  { 
    int k;
    for (k = 0; k < 2; k++)
      { r2_t *vk = &(edge->v[k]);
        fprintf(wr, "  v[%d] = ( %.6f %.6f )\n", k, vk->c[0], vk->c[1]);
      }
  }

i2_t stpoly_STP_round_point(r2_t *v, double eps, bool_t even)
  { 
    int rem = (even ? 0 : -1); /* Desired remainder, or {-1} is any. */
    i2_t qv;
    int k;
    for (k = 0; k < 2; k++)
      { qv.c[k] = (int32_t)iround(v->c[k], eps, 2, rem, INT32_MAX); }
    return qv;
  }

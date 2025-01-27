/* See argparser_extra.h. */
/* Last edited on 2025-01-22 19:06:27 by stolfi */

/* Copyright © 2003 Jorge Stolfi, Unicamp. See note at end of file. */
/* Based on Params.m3 by J.Stolfi, DEC-SRC, 1988.  */

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <float.h>

#include <affirm.h>
#include <jswsize.h>
#include <argparser.h>

#include <argparser_extra.h>

char *argparser_print_info_line(FILE *wr, char *lin, uint32_t wd);
  /* Prints the string {lin} (terminated by '\000' or '\n') to {wr}.
     Tries to breaks it so that it does not exceed {wd} characters per
     line, replicating its initial indentation on each subsequent
     piece. Such breaks are introduced at internal blanks. Any blanks
     at the break points or at the end of {lin} are discarded. Returns
     a pointer to the line's terminating character, which is is not
     printed. */
  
bool_t argparser_seems_keyword(char *n)
  { if (n[0] != '-') { return FALSE; }
    if (n[1] == '-') { return TRUE; }
    return (((n[1] >= 'a') && (n[1] <= 'z')) || ((n[1] >= 'A') && (n[1] <= 'Z')));
  }
  
bool_t argparser_key_matches(char *a, char *b)
  { if (strcmp(a, b) == 0) 
      { return TRUE; }
    else if ((a[0] == '-') && (b[0] == '-'))
      { if ((a[1] == '-') && (strcmp(a+1,b) == 0)) { return TRUE; }
        if ((b[1] == '-') && (strcmp(a,b+1) == 0)) { return TRUE; }
      }
    return FALSE;
  }
    
int64_t argparser_parse_int_string(argparser_t *pp, uint32_t index, char *arg, int64_t min, int64_t max, char **rest_P)
  { char *rest = NULL;
    errno = 0;
    int64_t v = (int64_t)strtoll(arg, &rest, 10);
    if (errno == EINVAL)
      { if (v == LONG_MAX)
          { argparser_error_at(pp, "integer overflow", "at", index); }
        else if (v == LONG_MIN)
          { argparser_error_at(pp, "integer underflow", "at", index); }
        else
          { assert(v == 0);
            argparser_error_at(pp, "invalid integer", "at", index);
          }
      }
    else if (rest == arg)
      { /* No digits parsed. Must be syntax error. */
        argparser_error_at(pp, "invalid integer", "at", index);
      }
    else if (((*rest) != '\000') && (rest_P == NULL))
      { /* No digits parsed. Must be syntax error. */
        argparser_error_at(pp, "invalid integer", "at", index);
      }
    else
      { if ((v > max) || (v < min)) 
         { argparser_arg_msg(pp, "", index, " should be in %s", "");
           fprintf(pp->wr, ("[%" int64_d_fmt "..%" int64_d_fmt "].\n"), min, max);
           argparser_print_help_and_halt(pp, 1);
         }
      }
    if (rest_P != NULL) { (*rest_P) = rest; }
    return v;
  }
 
void argparser_arg_msg(argparser_t *pp, char *msg1, uint32_t index, char *msg2, char *val)
  { fprintf(pp->wr, "%s: %s", pp->arg.e[0], msg1); 
    if ((index >= 1) && (index < pp->arg.ne))
      { fprintf(pp->wr, "parameter %d = \"%s\"", index, pp->arg.e[index]); }
    fprintf(pp->wr, msg2, val);
  }

void argparser_error_at(argparser_t *pp, char *msg, char* pos, uint32_t index)
  { fprintf(pp->wr, "%s: ** %s\n", pp->arg.e[0], msg);
    if ((index >= 1) && (index < pp->arg.ne))
      { fprintf(pp->wr, " %s argument %d = \"%s\"", pos, index, pp->arg.e[index]);
        fprintf(pp->wr, " (%s)", (pp->parsed.e[index] ? "parsed" : "unparsed"));
        fprintf(pp->wr, "\n");
      }
    argparser_print_help_and_halt(pp, 1);
  }

void argparser_print_info(argparser_t *pp, uint32_t wd)
  { for (uint32_t i = 0;  i < pp->nhelp; i++)
      { argparser_print_text(pp->wr, pp->info.e[i], wd); }
    fprintf(pp->wr, "\n");
  }

void argparser_print_text(FILE *wr, char *info, uint32_t wd)
  { char *lin = info;
    while ((*lin) != '\000')
      { /* Print line starting at {lin}, advance {lin} to its end: */
        lin = argparser_print_info_line(wr, lin, wd);
        /* Skip eol, if any: */
        if ((*lin) == '\n') { lin++; }
      }
  }

char *argparser_print_info_line(FILE *wr, char *lin, uint32_t wd)
  { /* Compute the indentation of this input line: */
    int32_t indent = 0;
    while ((*lin) == ' ') { lin++; indent++; }
    do 
      { /* Get next output line, from {lin} to {end}: */
        char *end = lin;
        char *q = end;
        while(((*q) != '\000') && ((*q) != '\n'))
          { /* Advance {q} to end of current word: */
            while (((*q) != ' ') && ((*q) != '\000') && ((*q) != '\n')) { q++; }
            /* Can we add that word to the output line? */
            if ((end > lin) && ((q - lin) + indent > wd)) { break; }
            /* Yes, add it: */
            end = q;
            /* Advance {q} to beginning of next word: */
            while ((*q) == ' ') { q++; }
          }
        /* Print the current output line with fixed indent: */
        fprintf(wr, "%*s%.*s\n", indent, "", (int32_t)(end-lin), lin);
        fflush(wr);
        lin = end;
        /* Skip blanks: */
        while ((*lin) == ' ') { lin++; }
      }
    while (((*lin) != '\000') && ((*lin) != '\n'));
    return lin;
  }

double_vec_t argparser_get_next_double_vec(argparser_t *pp, int32_t *NC)
  { double_vec_t v = double_vec_new(0);
    uint32_t NP = 0; /* Number of values actually parsed. */
    uint32_t NPMAX = ((NC == NULL) || ((*NC) < 0) ? UINT32_MAX : (uint32_t)(*NC)); /* Max to parse. */
    while ((NP < NPMAX) && argparser_next_is_number(pp))
      { double_vec_expand(&v, (vec_index_t)NP);
        v.e[NP] = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        NP++;
      }
    if ((NC != NULL) && (NP != 1))
      { /* Set/check {*NC} against {NP}: */
        if ((*NC) < 0)
          { (*NC) = (int32_t)NP; }
        else
          { if (NP != (uint32_t)(*NC)) { argparser_error(pp, "wrong number of elements"); } }
      }
    double_vec_trim(&v, NP);
    if (argparser_keyword_present_next(pp, "/"))
      { double den = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        if (den == 0.0) { argparser_error(pp, "bad denominator"); } 
        for (int32_t c = 0; c < NP; c++) { v.e[c] /= den; }
      }
    return v;
  }

int32_vec_t argparser_get_next_int32_vec(argparser_t *pp, int32_t *NC)
  { int32_vec_t v = int32_vec_new(0);
    uint32_t NP = 0; /* Number of values actually parsed. */
    uint32_t NPMAX = (uint32_t)((NC == NULL) || ((*NC) < 0) ? INT32_MAX : (*NC)); /* Max to parse. */
    while ((NP < NPMAX) && argparser_next_is_number(pp))
      { int32_vec_expand(&v, (vec_index_t)NP);
        v.e[NP] = (int32_t)argparser_get_next_int(pp, -INT32_MAX, +INT32_MAX);
        NP++;
      }
    if ((NC != NULL) && (NP != 1))
      { /* Set/check {*NC} against {NP}: */
        if ((*NC) < 0)
          { (*NC) = (int32_t)NP; }
        else
          { if (NP != (*NC)) { argparser_error(pp, "wrong number of elements"); } }
      }
    int32_vec_trim(&v, NP);
    return v;
  }
  
char *argparser_get_next_file_name(argparser_t *pp)
  { /* Peek at next argument: */
    char *nx = argparser_next(pp);
    /* If it is missing or parsed, it is not a file name: */
    if (nx == NULL) { return NULL; }
    /* If it is empty, it is not a file name: */
    if (nx[0] == '\000') { return NULL; }
    /* If it looks like a keyword, it is not a file name: */
    if ((nx[0] == '-') && (nx[1] != '\000')) { return NULL; }
    /* If it begins with a funny character, it is not a file name: */
    if
      ( (nx[0] != '@') &&
        (nx[0] != '/') && 
        (nx[0] != '.') &&
        ((nx[0] < 'A') || (nx[0] > 'Z')) &&
        ((nx[0] < 'a') || (nx[0] > 'z')) &&
        ((nx[0] < '0') || (nx[0] > '9'))
      ) 
      { return NULL; }
    /* Grab the next argument and check if it is a valid filename: */
    nx = argparser_get_next(pp);
    return nx;
  }

string_vec_t argparser_get_next_file_name_list(argparser_t *pp, uint32_t *NNP)
  { uint32_t NNMax = ((NNP == NULL) || ((*NNP) < 0) ? UINT32_MAX : (*NNP));
    uint32_t NN = 0;
    string_vec_t nvec = string_vec_new(0);
    while (NN < NNMax)
      { char *name = argparser_get_next_file_name(pp);
        if (name == NULL) { break; }
        string_vec_expand(&(nvec), (vec_index_t)NN);
        nvec.e[NN] = name;
        NN++;
      }
    if (NNP != NULL)
      { if ((*NNP) >= 0)
          { if (NN != (*NNP)) { argparser_error(pp, "not enough file names"); } }
        else
          { (*NNP) = NN; }
      }
    string_vec_trim(&(nvec), NN);
    return nvec;
  }

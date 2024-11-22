/* See fget_data.h */
/* Last edited on 2024-11-15 19:12:24 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <fget.h>
#include <bool.h>
#include <affirm.h>

#include <fget_data.h>

bool_t fget_data_fields(FILE *rd, char cmtc, int32_t nf, int8_t type[], char* alf[], double num[])
  {
    bool_t debug = FALSE;
    
    auto void debug_next(void);
    
    while(TRUE)
      { /* Skip spaces and comments on the current line: */
        if (debug) { fprintf(stderr, "checking for comment or endline...\n"); }
        bool_t eol = fget_test_comment_or_eol(rd, cmtc, NULL);
        if (eol) 
          { /* Consumed the end-of-line: */  
            if (debug) { fprintf(stderr, "skipped empty/blank/comment line\n"); }
            continue;
          }
        /* Stopped at non-comment, non-space before {eol}. Check it: */
        int32_t r = fgetc(rd); char c = (char)r;
        if (r == EOF) { return FALSE; }
        if (debug) 
          { fprintf(stderr, "found something that begins with '%c' = \\%03o\n", c, c);
            assert((c != '\000') && (c != ' ') && (c != '\240') && (c != '\011'));
          }
        assert(! fget_is_formatting_char(c));
        /* Not space, {cmtc}, eol, {EOF}: */
        ungetc(c, rd); 
        break;
      }
    /* Parse or skip the next {nf} data fields: */
    for (int32_t kf = 0; kf < nf; kf++)
      { alf[kf] = NULL; num[kf] = NAN;
        if (debug) { debug_next(); }
        if (type[kf] == 2)
          { /* Get a numeric field: */
            if (debug) { fprintf(stderr, "parsing field %d as numeric ...", kf); }
            /* This will skip spaces before the next field, gobble one */
            /* or more non-blank chars until the first char that is not */
            /* part of a number, then try to parse those chars as number: */
            num[kf] = fget_double(rd);
            if (debug) { fprintf(stderr, " got %24.16e\n", num[kf]); }
          }
        else
          { /* Get a non-numeric field or skip the next field: */
            if (debug) { fprintf(stderr, "parsing field %d as alpha ...", kf); }
            if (debug) { debug_next(); }
            fget_skip_spaces(rd);
            if (debug) { debug_next(); }
            /* This will gobble up zero of more chars until the */
            /* formatting char (including end-of-line) or {cmtc}): */
            char *fk = fget_to_delims(rd, cmtc, fget_formatting_chars);
            if (debug) { debug_next(); }
            /* If it hit a formatting char or {cmtc} right away, it*/
            /* would have returned an empty string: */
            assert(fk != NULL);
            if (debug) { fprintf(stderr, " got \"%s\"", fk); }
            demand ((*fk) != '\000', "data field not found");
            /* Store or ignore field: */
            if (type[kf] == 1)
              { if (debug) { fprintf(stderr, " (saved)\n"); }
                alf[kf] = fk;
              }
            else if (type[kf] == 0)
              { if (debug) { fprintf(stderr, " (discarded)\n"); }
                free(fk);
              }
            else
              { demand(FALSE, "invalid field type"); }    
          }
      }
    fget_skip_to_eol(rd);
      
    return TRUE;
        
    void debug_next(void)
      { 
        int32_t rr = fgetc(rd);
        if (rr == EOF) 
          { fprintf(stderr, " < at EOF >\n"); }
        else
          { fprintf(stderr, "  < next char = '\\%03o' = '%c' >\n", rr, (char)rr); ungetc(rr, rd); }
      }
  }

void fget_data_set_field_type(int32_t kf, int8_t tkf, bool_t rep_ok, int32_t nf, int8_t type[])
  { if (kf < 0) { return; }
    demand (kf < nf, "invalid data field index");
    demand((tkf >= 0) && (tkf <= 2), "invalid field type");
    if (type[kf] != 0)
      { if (tkf == type[kf])
          { demand(rep_ok, "data field cannot be used twice"); }
        else
          { demand(FALSE, "data field already declared as different type"); }
      }
    type[kf] = tkf;
  }

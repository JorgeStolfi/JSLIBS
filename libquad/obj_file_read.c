/* See {obj_file_read.h}. */
/* Last edited on 2024-12-22 10:46:32 by stolfi */
 
#define obj_file_read_C_copyright \
  "Copyright © 2024 State University of Campinas (UNICAMP).\n\n" jslibs_copyright
 
/* Written by J. Stolfi in June 2024. */ 

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <r3.h>
#include <affirm.h>
#include <fget.h>
#include <jsprintf.h>

#include <obj_file.h>
#include <obj_file_read.h>

#define debug FALSE

void obj_file_read_err(uint32_t nlin, char *msg);
  /* Prints the line number {nlin} and the message {msg}, then bombs out. */

void obj_file_read_command(FILE *rd, int32_t *ch1_P, int32_t *ch2_P, uint32_t *nlin_P);
  /* Assumes that the file {rd} is positioned at the beginning of a line.
    Skip spaces, then reads the next two characters, assumed to be the line type,
    and returns them in {*ch1_P} and {*ch2_P}.
    
    If the line type is a single letter ("v", "f", etc.) returns ' ' in {ch2}.
    
    If the file ends before it finds a line type, returns {EOF} in {ch1}.
    
    If the line type has more than two letters (like "usemtl" and "mtllib"),
    returns the first two letters in {ch1} and {ch2}, and skips the other letters.
    
    Failse if the face has fewer than 3 corners.
    
    Skips lines that are blank or contain only a "#"-comment. On input,
    {*nlin_P} should be the number of lines already parsed, not counting the 
    current one.  On output, {*nlin_P} is incremented by 1 plus the number of 
    blank/comment lines that were skipped. */

void obj_file_read_coords(FILE *rd, uint32_t nlin, char *elname, uint32_t ncmin, r3_t *u_P, char **lab_P);
  /* Reads coordinates of a vertex, texpoint, or normal {u} from {rd}, and returns 
    them as an {r3_t} triple in {*u_P}.  Fails if there are less than {ncmin} or more than 
    three coordinates on the line.  If fewer than 3 coordinates are read,
    the remaining ones are filled with zeros.
    
    If {lab_P} is not zero, returns in {*lab_P} a string that is the text of the '#'-comment
    on the data line that defined {u}; or a newly-allocated empty string {
    
    The parameter {nlin} should be the current line number in the file, and 
    {elname} should be "vertex", "normal", or "texpoint".  These
    parameters are used to form the error message if something goes wrong. */

int32_t obj_file_read_index(FILE *rd, uint32_t nlin, char *elname, uint32_t n); 
  /* If the next character in {rd} is space, '/', '#', or end-of-line, 
    leaves that character there and returns {-1}.  Otherwise tries to parse
    the next one or more characters as an integer {ix}, possibly with a minus sign.
    The procedure fails if it is not 
    
    If {ix} is positive, it must be between 1 and {n}, and the procedure
    returns {ix-1}.  if it is negative, it must be between {-n} and {-1},
    and the procedure returns {n+ix}.  Either way, the result will be in {0..n-1}.
    
    The parameter {nlin} should be the current line number in the file, and 
    {elname} should be "vertex", "normal", or "texpoint".  These
    parameters are used to form the error message if something goes wrong. */

void  obj_file_read_face
  ( FILE *rd,
    uint32_t nlin,
    uint32_t nv,
    int32_vec_t *FVk,
    uint32_t nt,
    int32_vec_t *FTk,
    uint32_t nn,
    int32_vec_t *FNk
  );
  /* Assumes that {rd} is positioned right after an "f" command. Reads the data for the corners
    of a face and stores them into the arrays {FVk}, {FTk}, {FNk}, expanding and trimming them as
    needed. 
    
    The rest of the line must be a sequence of corner specs.  A corner spec is the index
    of a vertex, a "/", the index of a texpoint, a "/", and the index of a normal.
    The last two maybe omitted, as well as any trailing "/" not followed by an index.
    
    The vertex index {ixv} must be in {1..nv}, or in {-nv..-1} indicating {nv+1-nv}.
    Likewise, the texpoint index {ixt} must be in {1..nt} or {-nt..-1},
    and the normal index must be in {1..nn} or {-nn..-1}.  Note that indices in the
    file start at 1 (not at 0, as they do in memory).
    
    On input the arrays {FVk,FTk,FNk} must be initialized with valid vectors, possibly empty.
    If they are not empty, the storage is reused and the old contents is overwritten.
    
    The parameter {nlin} should be the current input file 
    line number (from 1), for error messages. */

char *obj_file_read_cleanup_label(char *lab);
  /* If {lab} is not {NULL}, skips any leading '#' and spaces, as
    well as any trailing '.' and spaces. If the remaining string is
    not empty, returns a newly allocated copy of it. In all other
    cases returns {NULL}. */

obj_file_data_t *obj_file_read(FILE *rd, bool_t verbose)
  { uint32_t nlin = 0; /* Number of lines read. */
    uint32_t nv = 0; /* Number of vertices ('v' lines) found. */
    uint32_t nt = 0; /* Number of texpoints ('vt' lines) found. */
    uint32_t nn = 0; /* Number of normals ('vn' lines) found. */
    uint32_t nf = 0; /* Number of faces ('f' lines) found. */
    uint32_t ns = 0; /* Number of sides/corners (entries in 'f' lines) found. */
    
    obj_file_data_t *D = obj_file_data_new();
    
    while (TRUE)
      { int32_t ch1, ch2; /* The two characters of an OBJ line type. */
        obj_file_read_command(rd, &ch1, &ch2, &nlin);
        if (ch1 == EOF)
          { break; }
        if (debug) { fprintf(stderr, "    line %d command = \"%c%c\"\n", nlin, ch1, ch2); }
        if (ch1 == 'v')
          { /* Vertex, normal, or textmap line: */
            if (ch2 == ' ')
              { /* New vertex: */
                r3_t u = (r3_t){{ NAN, NAN, NAN }};
                char *lab = NULL;
                obj_file_read_coords(rd, nlin, "vertex", 3, &u, &lab);
                r3_vec_expand(&(D->V), (int32_t)nv);
                D->V.e[nv] = u;
                string_vec_expand(&(D->VL), (int32_t)nv);
                D->VL.e[nv] = obj_file_read_cleanup_label(lab);
                if (debug) { fget_show_next(stderr, "    next char = ", rd, "\n"); }
                fget_eol(rd);
                nv++;
              }
            else if (ch2 == 't')
              { /* New texpoint: */
                r3_t u = (r3_t){{ NAN, NAN, NAN }};
                obj_file_read_coords(rd, nlin, "texpoint", 2, &u, NULL);
                r3_vec_expand(&(D->T), (int32_t)nt);
                D->T.e[nt] = u;
                if (debug) { fget_show_next(stderr, "    next char = ", rd, "\n"); }
                fget_eol(rd);
                nt++;
              }
            else if (ch2 == 'n')
              { /* New normal vector: */
                r3_t u = (r3_t){{ NAN, NAN, NAN }};
                obj_file_read_coords(rd, nlin, "normal", 3, &u, NULL);
                double um = r3_dir(&u, &u); /* Normalize to unit length. */
                if (um < 1.0e-6) { obj_file_read_err(nlin, "normal vector is too short"); }
                r3_vec_expand(&(D->N), (int32_t)nn);
                D->N.e[nn] = u;
                if (debug) { fget_show_next(stderr, "    next char = ", rd, "\n"); }
                fget_eol(rd);
                nn++;
              }
            else
              { 
                char *msg = jsprintf("invalid line type \"%c%c\"", ch1, ch2);
                obj_file_read_err(nlin, msg);
              }
          }
        else if (ch1 == 'f')
          { /* Face definition line: */
            assert(ch2 == ' ');
            obj_file_face_vec_expand(&(D->FV), (int32_t)nf); D->FV.e[nf] = int32_vec_new(5); 
            obj_file_face_vec_expand(&(D->FT), (int32_t)nf); D->FT.e[nf] = int32_vec_new(5); 
            obj_file_face_vec_expand(&(D->FN), (int32_t)nf); D->FN.e[nf] = int32_vec_new(5); 
            obj_file_read_face(rd, nlin, nv, &(D->FV.e[nf]), nt, &(D->FT.e[nf]), nn, &(D->FN.e[nf]));
            ns += D->FV.e[nf].ne;
            nf++;
            if (debug) { fget_show_next(stderr, "    next char = ", rd, "\n"); }
            fget_eol(rd);
          }
        else if ((ch1 == 'o') && (ch2 == ' '))
          { /* Object name command. */
            fget_skip_to_eol(rd); 
          }
        else if ((ch1 == 'g') && (ch2 == ' '))
          { /* Grouping command: ignored. */
            fget_skip_to_eol(rd); 
          }
        else if ((ch1 == 'u') && (ch2 == 's'))
          { /* "usemtl" command: ignored. */
            fget_skip_to_eol(rd); 
          }
        else if ((ch1 == 'm') && (ch2 == 't'))
          { /* "mtllib" command: ignored. */
            fget_skip_to_eol(rd); 
          }
        else
          { char *msg = jsprintf("invalid line type \"%c%c\"", ch1, ch2);
            obj_file_read_err(nlin, msg);
          }
      }
    if (verbose) 
      { fprintf(stderr, "%s: found %d vertices  %d texpoints  %d normals  %d faces  %d corners\n", __FUNCTION__, nv, nt, nn, nf, ns); }
        
    r3_vec_trim(&(D->V), nv);                
    r3_vec_trim(&(D->T), nt);                
    r3_vec_trim(&(D->N), nn);                
    string_vec_trim(&(D->VL), nv);                
    obj_file_face_vec_trim(&(D->FV), nf);  
    obj_file_face_vec_trim(&(D->FT), nf);  
    obj_file_face_vec_trim(&(D->FN), nf); 
    
    return D;
  }

void obj_file_read_command(FILE *rd, int32_t *ch1_P, int32_t *ch2_P, uint32_t *nlin_P)
  { 
    /* Assumes that we are at beginning-of-line: */
    int32_t ch1 = '?', ch2 = '?';
    (*nlin_P)++;
    while (fget_test_comment_or_eol(rd, '#', NULL)) 
      { if (debug) { fprintf(stderr, "    line %d skipped\n", (*nlin_P)); }
        (*nlin_P)++;
      }
    if (fget_test_eof(rd)) 
      { ch1 = EOF; ch2 = ' '; }
    else
      { ch1 = fgetc(rd);
        assert((ch1 != '#') && (ch1 != '\n'));
        ch2 = fgetc(rd);
        if (ch2 == EOF) 
          { ch2 = ' '; }
        else if ((ch2 < 'a') || (ch2 > 'z')) 
          { ungetc(ch2, rd); ch2 = ' '; }
        else
          { assert((ch2 != '#') && (ch2 != '\n'));
            /* Skip rest of multi-letter command like "usemtl", "mtllib", "c_interp": */
            if ((ch1 >= 'a') && (ch1 <= 'z') && ((ch2 == '_') || ((ch2 >= 'a') && (ch2 <= 'z'))))
              { int32_t ch3;
                do { ch3 = fgetc(rd); } while ((ch3 == '_') || ((ch3 >= 'a') && (ch3 <= 'z')));
                if (ch3 != EOF) { ungetc(ch3, rd); }
              }
          }
      }
    (*ch1_P) = ch1;
    (*ch2_P) = ch2;
    return;
  }

void  obj_file_read_face
  ( FILE *rd,
    uint32_t nlin,
    uint32_t nv,
    int32_vec_t *FVk,
    uint32_t nt,
    int32_vec_t *FTk,
    uint32_t nn,
    int32_vec_t *FNk
  )  
  { uint32_t nc = 0;  /* Number of corners read for this face. */
    while(TRUE)
      { if (debug) { fprintf(stderr, "      looking for corner data, nc = %d\n", nc); }
        if (fget_test_comment_or_eol(rd, '#', NULL)) 
          { ungetc('\n', rd); break; }

        int32_t ixv = obj_file_read_index(rd, nlin, "vertex", nv);
        if (ixv == -1) { obj_file_read_err(nlin, "missing vertex index"); }

        int32_vec_expand(FVk, (int32_t)nc); 
        int32_vec_expand(FTk, (int32_t)nc); 
        int32_vec_expand(FNk, (int32_t)nc); 

        assert((ixv >= 1) && (ixv <= nv));
        FVk->e[nc] = ixv-1;
        FTk->e[nc] = -1;  /* If there is no '/'. */
        FNk->e[nc] = -1;  /* If there is no '/'. */

        if (fget_test_char(rd, '/'))
          { int32_t ixt = obj_file_read_index(rd, nlin, "texpoint", nt);
            if (ixt == -1)
              { FTk->e[nc] = -1; }
            else
              { assert((ixt >= 1) && (ixt <= nt));
                FTk->e[nc] = ixt-1;
              }
          }
        if (fget_test_char(rd, '/'))
          { int32_t ixn = obj_file_read_index(rd, nlin, "normal", nn);
            if (ixn == -1)
              { FNk->e[nc] = -1; }
            else
              { assert((ixn >= 1) && (ixn <= nn));
                FNk->e[nc] = ixn-1;
              }
          }
        nc++;
      }
    if (nc < 3)
      { 
        char *msg = jsprintf("face has only %d corners (min 3)", nc);
        obj_file_read_err(nlin, msg);
      }
    int32_vec_trim(FVk, nc);  
    int32_vec_trim(FTk, nc);  
    int32_vec_trim(FNk, nc);  
  }

int32_t obj_file_read_index(FILE *rd, uint32_t nlin, char *elname, uint32_t n)
  { 
    if (fget_test_char(rd, '/')) { ungetc('/', rd); return -1; }
    if (fget_test_char(rd, ' ')) { ungetc(' ', rd); return -1; }
    if (fget_test_comment_or_eol(rd, '#', NULL)) { ungetc('\n', rd); return -1; }
    int32_t ix = fget_int32(rd);
    if (ix < 0) 
      { /* Negative indices are relative to last item defined: */
        if (ix < -n) 
          { 
            char *msg = jsprintf("invalid negative %s index %d", elname, ix);
            obj_file_read_err(nlin, msg);
          }
        ix = (int32_t)n + 1 + ix;
      }
    else
      { if ((ix < 1) || (ix > n))
          { char *msg = jsprintf("invalid %s index %d", elname, ix);
            obj_file_read_err(nlin, msg);
          }
      }
    if (debug) { fprintf(stderr, "        read %s index %d\n", elname, ix); }
    if (debug) { int32_t chn = fgetc(rd); ungetc(chn, rd); fprintf(stderr, "      next char = '%c'\n", chn); }
    return ix;
  }

void obj_file_read_coords(FILE * rd, uint32_t nlin, char *elname, uint32_t ncmin, r3_t *u_P, char **lab_P)
  {
    r3_t u = (r3_t){{ 0,0,0 }};
    char *lab = NULL;
    uint32_t nc = 0; /* Number of coordinates read. */
    while (TRUE)
      { if (fget_test_comment_or_eol(rd, '#', &lab))
          { ungetc('\n', rd);
            if (debug) { fprintf(stderr, "      found '#' or LF, lab = \"%s\"\n", lab); }
            break;
          }
        double uj = fget_double(rd);
        if (debug) { fprintf(stderr, "      read %9.4f\n", uj); }
        if (debug) { fget_show_next(stderr, "    next char = ", rd, "\n"); }
        if (nc < 3) { u.c[nc] = uj; }
        nc++;
      }
    if (nc > 3) 
      { 
        char *msg = jsprintf("too many %s coordinates (%d, max 3)", elname, nc);
        obj_file_read_err(nlin, msg);
      }
    else if (nc < ncmin)
      { 
        char *msg = jsprintf("too few %s coordinates (%d, min %d)", elname, nc, ncmin);
        obj_file_read_err(nlin, msg);
      }
    (*u_P) = u;
    (*lab_P) = obj_file_read_cleanup_label(lab); free(lab);
    return;
  }
  
char *obj_file_read_cleanup_label(char *lab)
  { char *res = NULL;
    if (lab != NULL)
      { char *p = lab;
        while(((*p) != 0) && (((*p) == '#') || fget_is_space(*p))) { p++; }
        char (*q) = p + strlen(p);
        while((q > p) && ( fget_is_space(*q) || ((*q) == '.'))) { q--; }
        if (q > p) { 
          res = jsprintf("%*s", (int32_t)(q-p), p);
        }
      }
    return res;
  }

void obj_file_read_err(uint32_t nlin, char *msg)
  { fprintf(stderr, "line %d: ** %s: %s\n", nlin, __FUNCTION__, msg); 
    assert(FALSE);
  }

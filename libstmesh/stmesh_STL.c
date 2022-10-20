/* See {stmesh_STL.h} */
/* Last edited on 2022-10-20 06:03:10 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bool.h>
#include <fget.h>
#include <jsmath.h>
#include <bvtable.h>
#include <bvhash.h>
#include <argparser.h>

#include <stmesh_rep.h>
#include <stmesh_STL.h>

bool_t stmesh_STL_is_break(int32_t ch);
  /* TRUE iff {ch} is a line or page break character: ASCII CR ('\015'), 
    LF ('\012'), VT ('\013'), or FF ('\014'). */

bool_t stmesh_STL_is_space(int32_t ch);
  /* TRUE iff {ch} is a blank space character, but not a line or page
    break; namely, if {ch} is ASCII SP ('\040'), TAB ('\011'), NUL
    ('\000') and ISO Latin non-breaking space ('\240'). */
    
bool_t stmesh_STL_is_alpha(int32_t ch);
  /* TRUE iff {ch} is an upper or lower case ASCII letter, {'a'..'z'} or {'A'..'Z'}. */
    
bool_t stmesh_STL_is_digit(int32_t ch);
  /* TRUE iff {ch} is an ascii decimal digit, {'0'..'9'}. */
    
/* STL FILE PARSING

  The following procedures  increment {*lineP} whenever they skip over a line or page break character
  (except that Apple's CR-LF pair is counted as a single line break).  The {fileName}
  argument is used in error messages. */

bool_t stmesh_STL_skip_blanks(FILE *rd, char *fileName, int32_t *lineP);
  /* Skips space characters from {rd}, but not line or page breaks.  Returns TRUE if it finds
    a character other than those (which is not consumed), FALSE if it runs into
    the end-of-file.  */

bool_t stmesh_STL_skip_white_space(FILE *rd, char *fileName, int32_t *lineP);
  /* Skips space characters, line breaks, and page breaks from {rd}.  Returns TRUE if it finds
    a character other than those (which is not consumed), FALSE if it runs into
    the end-of-file.  */

char *stmesh_STL_get_keyword(FILE *rd, char *fileName, int32_t *lineP);
  /* Skips space characters, line breaks, and page breaks from {rd}. If
    it runs into end-of-file or a character that is not an ASCII letter,
    fails with an error message. Otherwise reads characters from {rd},
    until end-of-file, space, line break, or page break, and returns
    those characters as a newly allocated string. */

void stmesh_STL_check_keyword(FILE *rd, char *fileName, int32_t *lineP, char *key);
  /* Skips space characters, line breaks, and page breaks from {rd}
    and obtains the next token with {stmesh_STL_get_keyword}.
    Checks whether it is equal to {key}. Fails with an error message if
    there is no next token in {rd}, or the next token does not start
    with an ASCII letter, or it is not equal to {key}. */

char *stmesh_STL_parse_optional_name(FILE *rd, char *fileName, int32_t *lineP);
  /* Skips blanks (but not line breaks).  Then checks whether there is an 
    identifier in the current line; if there is, returns it as a freshly allocated string, and skips it;
    otherwise does nothing and returns {NULL}. */
    
bool_t stmesh_STL_ascii_read_face(FILE *rd, char *fileName, int32_t *lineP, stmesh_STL_face_t *face);
  /* Tries to read another triangular face from the ASCII STL file {rd}. If the
    next token in {rd} is "facet", reads the face data up to and
    including the "endfacet" token, stores that data in {*face}, and
    returns {TRUE}. If the next token is "endsolid" instead, consumes it
    and returns {FALSE}. Skips spaces, line breaks, and page breaks, and
    other formatting chars before and inside the parsed data. Fails with
    error message if the next token is something else, or the face
    is malformed. */

int32_t stmesh_STL_binary_read_header(FILE *rd, char *fileName, int32_t *lineP);
  /* Reads the 80-byte header of the binary STL file {rd}, and the 
    next line with the number of faces {nf}.  Returns {nf}.
    The line counter {*lineP} is incremented by 2. */

void stmesh_STL_binary_read_face(FILE *rd, char *fileName, int32_t *lineP, stmesh_STL_face_t *face);
  /* Tries to read another triangular face from the binary STL file {rd}. 
    Assumes that each face is 50 bytes: three components of the normal (3*float), 
    the coordinates of the three vertices (9*float), and a spacer (2 bytes).
    Increments the line counter {*lineP} by 1. */

/* IMPLEMENTATIONS */

bool_t stmesh_STL_is_break(int32_t ch)
  { 
    return (ch == '\012') || (ch == '\013') || (ch == '\014') || (ch == '\015');
  }

bool_t stmesh_STL_is_space(int32_t ch)
  { 
    return (ch == '\000') || (ch == '\011') || (ch == ' ') || (ch == '\240');
  }

bool_t stmesh_STL_is_alpha(int32_t ch)
  { 
    return ((ch >= 'a') && (ch <= 'z')) || ((ch >= 'A') && (ch <= 'Z'));
  }
    
bool_t stmesh_STL_is_digit(int32_t ch)
  { 
    return (ch >= '0') && (ch <= '9');
  }
    
bool_t stmesh_STL_skip_blanks(FILE *rd, char *fileName, int32_t *lineP)
  {
    int32_t ch = fgetc(rd);
    while (ch != EOF)
      { if (! stmesh_STL_is_space(ch))
          { ungetc(ch, rd); return TRUE; }
        ch = fgetc(rd);
      }
    return FALSE;
  }

bool_t stmesh_STL_skip_white_space(FILE *rd, char *fileName, int32_t *lineP)
  {
    int32_t prev_ch = EOF;  /* Previous character parsed in this call, or EOF if first. */
    int32_t ch = fgetc(rd);
    while (ch != EOF)
      { if (stmesh_STL_is_break(ch))
          { /* Possible line break, unless it is part of an Apple end-of-line (LF/VT/FF after a CR): */
            if ((prev_ch != '\015') || (ch != '\012')) { (*lineP)++; }
          }
        else if (! stmesh_STL_is_space(ch))
          { ungetc(ch, rd); return TRUE; }
        ch = fgetc(rd);
      }
    return FALSE;
  }
  
char *stmesh_STL_get_keyword(FILE *rd, char *fileName, int32_t *lineP)
  {
    if (! stmesh_STL_skip_white_space(rd, fileName, lineP))
      { fprintf(stderr, "%s:%d: ** expecting keyword, found end-of-file\n", fileName, (*lineP)); 
        exit(1);
      }
    int32_t ch = fgetc(rd);
    if (! stmesh_STL_is_alpha(ch))
      { fprintf(stderr, "%s:%d: ** expecting keyword, found '%c'\n", fileName, (*lineP), ch); 
        exit(1);
      }
    ungetc(ch, rd);
    return fget_string(rd);
  }

void stmesh_STL_check_keyword(FILE *rd, char *fileName, int32_t *lineP, char *key)
  {
    char *tok = stmesh_STL_get_keyword(rd, fileName, lineP);
    if (strcmp(tok, key) != 0)
      { fprintf(stderr, "%s:%d: ** expecting '%s', found '%s'\n", fileName, (*lineP), key, tok);
        exit(1);
      }
    free(tok);
  }

char *stmesh_STL_parse_optional_name(FILE *rd, char *fileName, int32_t *lineP)
  {
    if (! stmesh_STL_skip_blanks(rd, fileName, lineP))
      { return NULL; }
    int32_t ch = fgetc(rd);
    ungetc(ch, rd);
    if (stmesh_STL_is_alpha(ch))
      { char *s = fget_string(rd);
        return s;
      }
    else
      { return NULL; }
  }

void stmesh_STL_gen_read(char *fileName, bool_t binary, stmesh_STL_face_proc_t *process_face)  
  {
    char *ftype = (binary ? "rb" : "r");
    FILE *rd = fopen(fileName, ftype);
    if (rd == NULL) 
      { fprintf(stderr, "** failed to open file '%s'\n", fileName);
        exit(1);
      }

    stmesh_STL_face_t face;
    int32_t line = 1; /* Line number in file. */
    int32_t nf = 0;   /* Number of faces read from the STL file. */
    
    if (binary)
      { int32_t nf = stmesh_STL_binary_read_header(rd, fileName, &line);
        int32_t it;
        for (it = 0; it < nf; it++) 
          { stmesh_STL_binary_read_face(rd, fileName, &line, &face);
            process_face(line, &face);
          }
      }
    else
      { stmesh_STL_check_keyword(rd, fileName, &line, "solid");
        /* Some files have a solid name following the keyword: */
        char *name = stmesh_STL_parse_optional_name(rd, fileName, &line);
        if (name != NULL) { free(name); }
        /* Read the faces. */
        while (stmesh_STL_ascii_read_face(rd,fileName, &line, &face))
          { nf++;
            process_face(line, &face);
          }
      }

    fclose(rd);
  }

bool_t stmesh_STL_ascii_read_face(FILE *rd, char *fileName, int32_t *lineP, stmesh_STL_face_t *face)
  { char *tok = stmesh_STL_get_keyword(rd, fileName, lineP);
    if (strcmp(tok, "endsolid") == 0)
      { return FALSE; }
    else if (strcmp(tok, "facet") == 0)
      { int32_t i, k;
        stmesh_STL_check_keyword(rd, fileName, lineP, "normal");
        for (i = 0; i < 3; i++) 
          { if (! stmesh_STL_skip_white_space(rd, fileName, lineP))
              { fprintf(stderr, "%s:%d: ** expecting normal component, found end-of-file\n", fileName, (*lineP));
                exit(1);
              }
            face->normal.c[i] = (float)fget_double(rd);
          }
        stmesh_STL_check_keyword(rd, fileName, lineP, "outer"); 
        stmesh_STL_check_keyword(rd, fileName, lineP, "loop"); 
        for (k = 0; k < 3; k++) 
          { stmesh_STL_check_keyword(rd, fileName, lineP, "vertex"); 
            for (i = 0; i < 3; i++) 
              { if (! stmesh_STL_skip_white_space(rd, fileName, lineP))
                  { fprintf(stderr, "%s:%d: ** expecting vertex coord, found end-of-file\n", fileName, (*lineP));
                    exit(1);
                  }
                face->v[k].c[i] = (float)fget_double(rd);
              }
          }
        stmesh_STL_check_keyword(rd, fileName, lineP, "endloop"); 
        stmesh_STL_check_keyword(rd, fileName, lineP, "endfacet");
        return TRUE;
      }
    else
      { fprintf(stderr, "%s:%d: ** expected 'facet' or 'endsolid', found '%s'\n", fileName, (*lineP), tok);
        exit(1);
      }
  }

int32_t stmesh_STL_binary_read_header(FILE *rd, char *fileName, int32_t *lineP)
  { 
    char title[80];
    int32_t nf;
    int32_t err;
    err = (int32_t)fread(title, 80, 1, rd);
    if (err != 1) { fprintf(stderr, "%s: error reading binary STL file header\n", fileName); exit(1); }
    (*lineP)++;
    err = (int32_t)fread((void*)(&nf), 4, 1, rd);
    if (err != 1) { fprintf(stderr, "%s: error reading {nf} from binary STL file\n", fileName); exit(1); }
    (*lineP)++;
    return nf;
  }

void stmesh_STL_binary_read_face(FILE *rd, char *fileName, int32_t *lineP, stmesh_STL_face_t *face)
  { 
    /* Read the coordinates of normal and 3 vertices: */
    float vc[12];  /* Coordinates of normal (3) and vertices (3*3). */
    int32_t err;
    int32_t k;
    for (k = 0; k < 12; k++) 
      { err = (int32_t)fread((void*)(&vc[k]), sizeof(float), 1, rd);
        if (err != 1) { fprintf(stderr, "%s:%d: error reading binary STL datum %d\n", fileName, (*lineP), k); exit(1); }
      }
    
    /* Read the 16-bit padding: */
    unsigned short uint16;  /* Padding between faces. */
    err = (int32_t)fread((void*)(&uint16), sizeof(unsigned short), 1, rd); // spacer between successive faces
    if (err != 1) { fprintf(stderr, "%s:%d: error reading binary STL padding\n", fileName, (*lineP)); exit(1); }
    (*lineP)++;
    
    /* Fill the {stmesh_STL_face_t}. */
    face->normal = (stmesh_STL_r3_t){{  vc[0],  vc[1],  vc[2] }};
    face->v[0] = (stmesh_STL_r3_t){{  vc[3],  vc[4],  vc[5] }};
    face->v[1] = (stmesh_STL_r3_t){{  vc[6],  vc[7],  vc[8] }};
    face->v[2] = (stmesh_STL_r3_t){{  vc[9], vc[10], vc[11] }};
  }

void stmesh_STL_print_triangle(FILE *wr, stmesh_STL_face_t *f)
  { 
    int32_t k;
    for (k = 0; k < 3; k++)
      { stmesh_STL_r3_t *vk = &(f->v[k]);
        fprintf(wr, "  v[%d] = ( %.8f %.8f %.8f )\n", k, vk->c[0], vk->c[1], vk->c[2]);
      }
    stmesh_STL_r3_t *vn = &(f->normal);
    fprintf(wr, "  norm = ( %.8f %.8f %.8f )\n", vn->c[0], vn->c[1], vn->c[2]);
  }
 
i3_t stmesh_STL_round_point(stmesh_STL_r3_t *v, float eps, bool_t even)
  { 
    int32_t rem = (even ? 0 : -1); /* Desired remainder, or {-1} is any. */
    i3_t qv;
    int32_t k;
    for (k = 0; k < 3; k++)
      { qv.c[k] = (int32_t)iroundfrac(v->c[k], eps, 2, rem, INT32_MAX); }
    return qv;
  }

stmesh_t stmesh_STL_read(char *fileName, bool_t binary, float eps, uint32_t nfGuess, bool_t even, bool_t checkSorted)
  {
    char *format = ((char *[2]){ "ascii", "binary" })[binary];
    fprintf(stderr, "reading mesh from file %s (%s)\n", fileName, format);
    fprintf(stderr, "expecting about %u triangles\n", nfGuess);
    fprintf(stderr, "quantizing vertex coords to%s multiples of %.8f mm\n", (even ? " even" : ""), eps);
    if (checkSorted) { fprintf(stderr, "expecting triangles sorted by {.minZ}\n"); }

    /* While reading the STL file and building the structure, the
      entries must be linked by *indices* of elements, rather
      than pointers, because the hash tables may be reallocated. Only
      after the reading is finished can we translate the indices into
      pointers.
      
      During the reading, a vertex is represented by its
      integer position (as an {i3.t}); an (unoriented) edge is represented
      by indices of the two endpoints, in increasing order (as an
      {stmesh_vert_unx_pair_t}); and a face is represented by the indices of the three
      (unoriented) edges that are its sides (as an {stmesh_edge_unx_triple_t}),
      also in increasing order.
      
      The fields {uf.minZ} and {uf.maxZ}, that are not needed during reading,
      are set later, when these pairs and triples of integers are converted to
      records of types {stmesh_vert_rep_t}, {stmesh_edge_rep_t}, and
      {stmesh_face_rep_t}. */

    /* Lookup table to uniquify faces: */
    size_t sz_face = sizeof(stmesh_edge_unx_triple_t);
    uint32_t ng_face = nfGuess;  /* Expected number of faces. */
    bvtable_t *tb_face = bvtable_new(sz_face, ng_face);
    auto uint64_t hash_face(void *fp, size_t sz);
    auto int32_t cmp_face(void *xp, void *yp, size_t sz);
    
    /* Lookup table to uniquify edges: */
    size_t sz_edge = sizeof(stmesh_vert_unx_pair_t);
    uint32_t ng_edge = (3 * ng_face + 2)/2;
    bvtable_t *tb_edge = bvtable_new(sz_edge, ng_edge);
    auto uint64_t hash_edge(void *ep, size_t sz);
    auto int32_t cmp_edge(void *xp, void *yp, size_t sz);

    /* Lookup table to uniquify vertices: */
    size_t sz_vert = sizeof(i3_t);
    uint32_t ng_vert = (ng_face + 1)/2;
    bvtable_t *tb_vert = bvtable_new(sz_vert, ng_vert);
    auto uint64_t hash_vert(void *vp, size_t sz);
    auto int32_t cmp_vert(void *xp, void *yp, size_t sz);
    
    auto void process_STL_face(int32_t line, stmesh_STL_face_t *face);
      /* Procedure that quantizes an STL face {face} and stores it in
         the mesh, if not degenerate. */

    uint32_t nf_read = 0;  /* Number of faces read from the STL file. */
    uint32_t nf_keep = 0;  /* Number of faces retained in the mesh. */
    
    int32_t prevZ = INT32_MIN; /* If {checkSorted}, the {minZ} of the next face must be at least this big. */

    stmesh_STL_gen_read(fileName, binary, &process_STL_face);
    
    fprintf(stderr, "read %u triangles, kept %u\n", nf_read, nf_keep);

    /* Close the tables and get the basic data: */
    uint32_t nv, ne, nf; /* Number of vertices, edges, and faces. */
    i3_t *vpos;                     /* Quantized coordinates of vertices. */
    stmesh_vert_unx_pair_t *endv;   /* Edges, as pairs of vertex indices. */
    stmesh_edge_unx_triple_t *side; /* Faces, as triplets of unoriented edge indices. */
    
    bvtable_close(tb_vert, &nv, (void**)&(vpos)); 
    bvtable_close(tb_edge, &ne, (void**)&(endv)); 
    bvtable_close(tb_face, &nf, (void**)&(side)); 
    
    assert(nf == nf_keep);
    fprintf(stderr, "found %u distinct vertices and %u distinct edges\n", nv, ne);

    /* Build the mesh data structure. */
    stmesh_t mesh = stmesh_build(eps, nv, vpos, ne, endv, nf, side, checkSorted);
    
    return mesh;
    
    /* INTERNAL IMPLEMENTATIONS */
    
    void process_STL_face(int32_t line, stmesh_STL_face_t *stl_face)
      { 
        nf_read++;
        
        /* Quantize vertices, assign indices: */
        stmesh_edge_unx_triple_t uxv; /* Indices of corner vertices, assigned or recovered. */
        int32_t minZ = INT32_MAX; /* Minimum {Z}-coordinate, to check order. */
        int32_t k;
        for (k = 0; k < 3; k++)
          { /* Quantize the coordinates of corner {k}: */
            i3_t vposk = stmesh_STL_round_point(&(stl_face->v[k]), eps, even);
            /* Update the min {Z} coordinate of the triangle: */
            int32_t zk = vposk.c[2];
            if (zk < minZ) { minZ = zk; }
            /* Get the unique vertex index: */
            uxv.c[k] = bvtable_add(tb_vert, (void*)(&vposk), &hash_vert, &cmp_vert);
            demand(uxv.c[k] <= stmesh_nv_MAX, "too many vertices in mesh");
            /* Check for repeated vertices: */
            int32_t i;
            for (i = 0; i < k; i++)
              { if (uxv.c[i] == uxv.c[k])
                  { fprintf(stderr, "%s:%d: !! warning: vertices %d %d coincide, triangle ignored\n", fileName, line, i, k);
                    stmesh_STL_print_triangle(stderr, stl_face);
                    fprintf(stderr, "\n");
                    return;
                  }
              }
          }

        if (checkSorted)
          { /* Check that the {.minZ} fields are non-decreasing: */
            if (minZ < prevZ)
              { fprintf(stderr, "%s:%d: ** error: triangles not sorted by {minZ}\n", fileName, line);
                exit(1);
              }
            prevZ = minZ;
          }
        
        /* Get or assign indices to edges: */
        stmesh_edge_unx_triple_t uxside; /* Indices of the three unoriented edges that bound the face. */
        for (k = 0; k < 3; k++)
          { stmesh_vert_unx_pair_t uxendvk;       /* Temporary record for edge {k}. */
            stmesh_vert_unx_t uxv0 = uxv.c[k];       /* One endpoint. */
            stmesh_vert_unx_t uxv1 = uxv.c[(k+1)%3]; /* The other endpoint. */
            assert(uxv0 != uxv1);         /* We should have discarded degenerate triangles before. */
            /* Make sure that {uxv0} is the vertex with lowest index: */
            if (uxv0 > uxv1) { stmesh_vert_unx_t t = uxv0; uxv0 = uxv1; uxv1 = t; }
            uxendvk.c[0] = uxv0;
            uxendvk.c[1] = uxv1;
            uxside.c[k] = bvtable_add(tb_edge, (void*)(&uxendvk), &hash_edge, &cmp_edge);
            demand(uxside.c[k] <= stmesh_ne_MAX, "too many edges in mesh");
          }
          
        /* Sort the edge indices in increasing order: */
        if (uxside.c[0] > uxside.c[1]) { stmesh_edge_unx_t t = uxside.c[0]; uxside.c[0] = uxside.c[1]; uxside.c[1] = t; }
        if (uxside.c[1] > uxside.c[2]) { stmesh_edge_unx_t t = uxside.c[1]; uxside.c[1] = uxside.c[2]; uxside.c[2] = t; }
        if (uxside.c[0] > uxside.c[1]) { stmesh_edge_unx_t t = uxside.c[0]; uxside.c[0] = uxside.c[1]; uxside.c[1] = t; }
        
        /* Since the three vertices are distinct, the three sides must be distinct edges too. */

        /* Assign an index to the face: */
        stmesh_face_unx_t uxf = bvtable_add(tb_face, (void *)(&uxside), &hash_face, &cmp_face);
        if (uxf < nf_keep) 
          { /* Repeated face: */
            fprintf(stderr, "%s:%d: !! repeated triangle, ignored\n", fileName, line);
          }
        else
          { nf_keep++; }
      }

    /* Hashing and comparison procedures for faces, edges, and vertices: */
    
    uint64_t hash_face(void *fp, size_t sz)
      { assert(sz == sizeof(stmesh_edge_unx_triple_t));
        stmesh_edge_unx_triple_t *f = (stmesh_edge_unx_triple_t *)fp;
        /* Requires the edge indices {f.c[0..2]} to be increasing: */
        assert((f->c[0] < f->c[1]) && (f->c[1] < f->c[2]));
        /* Hash side indices: */
        uint64_t h = bvhash_bytes(fp, sizeof(stmesh_edge_unx_triple_t));
        return h;
      }
        
    auto int32_t cmp_face(void *xp, void *yp, size_t sz)
      { assert(sz == sizeof(stmesh_edge_unx_triple_t));
        stmesh_edge_unx_triple_t *x = (stmesh_edge_unx_triple_t *)xp;
        stmesh_edge_unx_triple_t *y = (stmesh_edge_unx_triple_t *)yp;
        /* Compare the edge indices lexicographically: */
        int32_t k;
        for (k = 0; k < 3; k++)
          { stmesh_edge_unx_t uxex = x->c[k];
            stmesh_edge_unx_t uxey = y->c[k];
            if (uxex < uxey)
              { return -1; }
            else if (uxex > uxey)
              { return +1; }
          }
        return 0;
      }
        
    uint64_t hash_edge(void *ep, size_t sz)
      { assert(sz == sizeof(stmesh_vert_unx_pair_t));
        stmesh_vert_unx_pair_t *e = (stmesh_vert_unx_pair_t *)ep;
        /* Requires endpoint indices in increasing order: */
        assert(e->c[0] < e->c[1]);
        /* Hash the endpoint indices: */
        uint64_t h = bvhash_bytes(ep, sizeof(stmesh_vert_unx_pair_t));
        return h;
      }
        
    auto int32_t cmp_edge(void *xp, void *yp, size_t sz)
      { assert(sz == sizeof(stmesh_vert_unx_pair_t));
        stmesh_vert_unx_pair_t *x = (stmesh_vert_unx_pair_t *)xp;
        stmesh_vert_unx_pair_t *y = (stmesh_vert_unx_pair_t *)yp;
        /* Compare the endpoint indices lexicographically: */
        int32_t k;
        for (k = 0; k < 2; k++)
          { stmesh_vert_unx_t uxvx = x->c[k];
            stmesh_vert_unx_t uxvy = y->c[k];
            if (uxvx < uxvy)
              { return -1; }
            else if (uxvx > uxvy)
              { return +1; }
          }
        return 0;
      }
        
    uint64_t hash_vert(void *vp, size_t sz)
      { assert(sz == sizeof(i3_t));
        /* Hash the quantized coords: */
        uint64_t h = bvhash_bytes(vp, sizeof(i3_t));
        return h;
      }
        
    auto int32_t cmp_vert(void *xp, void *yp, size_t sz)
      { assert(sz == sizeof(i3_t));
        i3_t *x = (i3_t *)xp;
        i3_t *y = (i3_t *)yp;
        /* Compare quantized coords lexicographically in order {Z,Y,X}: */
        int32_t k;
        for (k = 0; k < 3; k++)
          { int32_t xk = x->c[2-k];
            int32_t yk = y->c[2-k];
            if (xk < yk)
              { return -1; }
            else if (xk > yk)
              { return +1; }
          }
        return 0;
      }

  }

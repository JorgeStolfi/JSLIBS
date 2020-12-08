/* See {stmesh_section.h}. */
/* Last edited on 2015-11-16 03:00:31 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <affirm.h>
#include <r2.h>
#include <vec.h>

#include <stmesh.h>
#include <stmesh_section.h>

/* INTERNAL PROTOTYPES */

void stmesh_section_write_header(FILE *wr, stmesh_section_t *sec);
  /* Writes to {wr} the header and preamble of a file that will
    contain the cross-section {sec}. */
      
void stmesh_section_write_footer(FILE *wr);
  /* Writes to {wr} the footer of a file that contains a mesh cross-section. */
 
void stmesh_section_write_path(FILE *wr, r2_t v[], uint32_t ini, uint32_t fin);
  /* Writes to {wr} an open or closed path that is part of a cross-section.
    The vertices of the path are {v[ini..fin]}. */


void stmesh_section_read_header(FILE *rd, stmesh_section_t *sec);
  /* Reads from {rd} the header and preamble of a cross-section
    file, returning the fundamental length in {sec.eps}, the quantized {Z}-coordinate
    of the slicing plane in {sec.pZ}, the total number of segments in {sec.ns},
    and the number of paths {sec->nc}.  The tables
    {sec->v} and {sec->cstart} are set to NULL.
    
    Skips blank lines before the header line "begin stmesh_section". */

void stmesh_section_read_footer(FILE *drdr);
  /* Parses the footer line "end stmesh_section", including the final newline. */ 

void stmesh_section_read_path(FILE *rd, r2_t v[], uint32_t ini, uint32_t nv_max, uint32_t *finP);
  /* Reads another open or closed path from section file {rd}, from
    the header line "begin path" to the footer line "end path".
    Skips blank lines before the header, and consumes the end-of-line
    after the footer.
    
    The points are returned in {v[ini..fin]} where {ini} is the specified
    index and {fin} is determined by the procedure. The number
    of vertices in the path {fin+1-ini} will be the {nPoints} field
    after "begin path".  The procedure assumes that {v} was
    allocated with {nv_max} elements.  */


/* IMPLEMENTATIONS */

stmesh_section_t* stmesh_section_make(stmesh_t mesh, int32_t pZ, uint32_t nc, uint32_t estart[], stmesh_edge_t e[])
  { 
    bool_t debug = FALSE;
    
    demand(estart[0] == 0, "bad {estart}");
    uint32_t nv = estart[nc];
    demand(nv >= nc, "bad {nv,nc}");

    stmesh_section_t *sec = notnull(malloc(sizeof(stmesh_section_t)), "no mem");
    
    sec->eps = stmesh_get_eps(mesh);
    sec->pZ = pZ;
    sec->ns = nv - nc;
    sec->nc = nc;
    sec->v = notnull(malloc(nv*sizeof(r2_t)), "no mem");
    sec->cstart = notnull(malloc((nc+1)*sizeof(uint32_t)), "no mem");
    
    if (debug) { fprintf(stderr, "section with %d segments in %d paths with %d vertices\n", sec->ns, sec->nc, nv); }
    
    uint32_t iv = 0;
    sec->cstart[0] = iv;
    uint32_t kc;
    for (kc = 0; kc < nc; kc++)
      { 
        demand(estart[kc] == iv, "inconsistent {estart}");
        uint32_t fin = estart[kc+1] - 1;
        demand(fin >= iv+1, "empty or trivial path");
        if (debug) { fprintf(stderr, "  path %d vertices %d..%d\n", kc, iv, fin); }
        while (iv <= fin)
          { r2_t vi = stmesh_edge_plane_intersection(e[iv], pZ, sec->eps);
            sec->v[iv] = vi;
            iv++;
          }
        sec->cstart[kc+1] = iv;
      }
    if (debug) { fprintf(stderr, "  total %d vertices\n", iv); }
    assert(iv == nv);
    return sec;
  }

void stmesh_section_free(stmesh_section_t *sec)
  { 
    free(sec->v);
    free(sec->cstart);
    free(sec);
  }

void stmesh_section_write(FILE *wr, stmesh_section_t *sec)
  { 
    stmesh_section_write_header(wr, sec);
    
    assert(sec->cstart[0] == 0);
    uint32_t ini = 0; /* Start of next path. */
    uint32_t kc;
    for (kc = 0; kc < sec->nc; kc++)
      { uint32_t fin = sec->cstart[kc+1] - 1;
        assert(fin >= ini+1);
        stmesh_section_write_path(wr, sec->v, ini, fin);
        ini = fin + 1;
      }
    assert(ini == sec->ns + sec->nc);

    stmesh_section_write_footer(wr);
  }

void stmesh_section_write_header(FILE *wr, stmesh_section_t *sec)
  { 
    float eps = sec->eps;
    filefmt_write_header(wr, "stmesh_section", stmesh_section_VERSION);
    fprintf(wr, "eps = %23.15e\n", (double)eps);
    fprintf(wr, "planeZ = %+23.15e\n", ((double)sec->pZ)*(double)eps);
    fprintf(wr, "planeQ = %+11d\n", sec->pZ);
    fprintf(wr, "nSegs = %u\n", sec->ns);
    fprintf(wr, "nPaths = %u\n", sec->nc);
    fflush(wr);
  }

void stmesh_section_write_footer(FILE *wr)
  { 
    filefmt_write_footer(wr, "stmesh_section");
  }

void stmesh_section_write_path(FILE *wr, r2_t v[], uint32_t ini, uint32_t fin)
  { 
    uint32_t nv = fin - ini + 1; /* Number of points in path. */
    
    fprintf(wr, "\n");
    fprintf(wr, "begin path\n");
    fprintf(wr, "nPoints = %u\n", nv);
    
    uint32_t iv;
    for (iv = ini; iv <= fin; iv++)
      { r2_t *vi = &(v[iv]);
        fprintf(wr, " %+.3f %+.3f\n", (float)vi->c[0], (float)vi->c[1]); 
      }
    
    fprintf(wr, "end path\n");
  }

stmesh_section_t* stmesh_section_read(FILE *rd)
  {
    stmesh_section_t *sec = notnull(malloc(sizeof(stmesh_section_t)), "no mem");
    stmesh_section_read_header(rd, sec);
    
    /* Allocate vectors, to be trimmed later: */
    uint32_t nv_max = 2*sec->ns; /* Upper bound to the number of vertices in all paths. */
    uint32_t nc_max = sec->ns;   /* Upper bound to the number of paths. */
    sec->v = notnull(malloc(nv_max*sizeof(r2_t)), "no mem");
    sec->cstart = notnull(malloc((nc_max+1)*sizeof(uint32_t)), "no mem");

    sec->cstart[0] = 0;
    uint32_t ini = 0; /* Start of next path. */
    uint32_t kc;
    for (kc = 0; kc < sec->nc; kc++)
      { uint32_t fin;
        stmesh_section_read_path(rd, sec->v, ini, nv_max, &fin);
        assert(fin >= ini+1);
        ini = fin + 1;
        sec->cstart[kc+1] = ini;
      }
    assert(ini == sec->ns + sec->nc);

    stmesh_section_read_footer(rd);
    
    return sec;
  }

void stmesh_section_read_header(FILE *rd, stmesh_section_t *sec)
  { 
    /* Read the first {filefmt} header line and check the version: */
    filefmt_read_header(rd, "stmesh_section", stmesh_section_VERSION);
    
    /* Parse the preamble fields: */
    float eps = (float)nget_double(rd, "eps"); fget_eol(rd);
    double fpZ = nget_double(rd, "planeZ"); fget_eol(rd);
    int32_t pZ = nget_int32(rd, "planeQ"); fget_eol(rd);
    uint32_t nSegs = nget_uint32(rd, "nSegs", 10); fget_eol(rd);
    uint32_t nPaths = nget_uint32(rd, "nPaths", 10); fget_eol(rd);
    demand(eps > 0.0, "the {eps} field must be positive");
    demand(fabs(fpZ/eps - (double)pZ) < 0.0001, "inconsitent {planeQ} in file");
    demand(nPaths <= nSegs, "invalid num of paths");
    
    sec->eps = eps;
    sec->pZ = pZ;
    sec->ns = nSegs;
    sec->nc = nPaths;
    sec->v = NULL;
    sec->cstart = NULL;
  }

void stmesh_section_read_footer(FILE *rd)
  { 
    filefmt_read_footer(rd, "stmesh_section");
  }

void stmesh_section_read_path(FILE *rd, r2_t v[], uint32_t ini, uint32_t nv_max, uint32_t *finP)
  { 
    /* Read the line "begin path": */
    filefmt_read_header(rd, "path", NULL);

    /* Get the number of points of this path/loop: */
    uint32_t nPoints = nget_uint32(rd, "nPoints", 10); fget_eol(rd);
    demand(nPoints >= 2, "invalid number of points in path");
    uint32_t fin = ini + nPoints - 1;
    demand(fin < nv_max, "too many points in section");

    fget_skip_formatting_chars(rd);
    
    uint32_t iv;
    for(iv = ini; iv <= fin; iv++)
      { r2_t *vi = &(v[iv]);
        vi->c[0] = fget_double(rd);
        vi->c[1] = fget_double(rd);
        fget_eol(rd);
      }
    (*finP) = fin;
      
    /* Read the line "end path": */
    filefmt_read_footer(rd, "path");
  }
      
        
        

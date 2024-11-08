/* See argparser_geo.h. */
/* Last edited on 2024-10-31 02:46:10 by stolfi */

/* Copyright © 2003 Jorge Stolfi, Unicamp. See note at end of file. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <argparser.h>

#include <argparser_geo.h>

void argparser_get_next_rn(argparser_t *pp, double p[], int32_t n, double min, double max)
  { int32_t i;
    for (i = 0; i < n; i++) { p[i] = argparser_get_next_double(pp, min, max); }
  }

r2_t argparser_get_next_r2(argparser_t *pp, double min, double max)
  { r2_t p;
    argparser_get_next_rn(pp, p.c, 2, min, max);
    return p;
  } 

r3_t argparser_get_next_r3(argparser_t *pp, double min, double max)
  { r3_t p;
    argparser_get_next_rn(pp, p.c, 3, min, max);
    return p;
  } 

r4_t argparser_get_next_r4(argparser_t *pp, double min, double max)
  { r4_t p;
    argparser_get_next_rn(pp, p.c, 4, min, max);
    return p;
  } 
  
r6_t argparser_get_next_r6(argparser_t *pp, double min, double max)
  { r6_t p;
    argparser_get_next_rn(pp, p.c, 6, min, max);
    return p;
  } 
  
r3_t argparser_get_next_r3_dir(argparser_t *pp)
  { r3_t d;
    argparser_get_next_rn(pp, d.c, 3, -DBL_MAX, +DBL_MAX);
    r3_dir(&d, &d);
    return d;
  } 
 
void argparser_get_next_adjust(argparser_t *pp, double *adjP, double min, double max)
  { if (adjP != NULL) 
      { if (argparser_keyword_present_next(pp, "adjust"))
          { (*adjP) = argparser_get_next_double(pp, min, max); }
        else
          { (*adjP) = 0; }
      }
  }

hr2_pmap_type_t argparser_get_next_hr2_pmap_type(argparser_t *pp)
  { char *tname = argparser_get_next_non_keyword(pp);
    hr2_pmap_type_t type;
    bool_t ok = hr2_pmap_type_from_string(tname, &type);
    if (ok)
      { return type; }
    else
      { argparser_error(pp, "invalid projective map type");
        assert(FALSE); /* Should not get here. */
      }
  }

#define BIG_ELEM (1.0e+100)
  /* A very large matrix element value, but still far from overflow. */

hr2_pmap_t argparser_get_next_proj_map_matrix(argparser_t *pp)
  { /* Parse the nine entries of the projective matrix, row order: */
    r3x3_t A;
    int32_t i, j;
    for (i = 0; i < 3; i++)
      { for (j = 0; j < 3; j++)
          { A.c[i][j] = argparser_get_next_double(pp, -BIG_ELEM, +BIG_ELEM); }
      }
    /* Assemble the projective map: */
    hr2_pmap_t M;
    M.dir = A;
    r3x3_inv(&A, &(M.inv));
    return M;
  }

#define BIG_COORD (1.0e+100)
  /* A very large coordinate value, but still far from overflow. */

hr2_pmap_t argparser_get_next_proj_map_from_points(argparser_t *pp)
  { /* Parse the four input points {ip[0..3]}: */
    hr2_point_t ip[4];
    int32_t i;
    for (i = 0; i < 4; i++)
      { ip[i].c.c[0] = 1.0;
        ip[i].c.c[1] = argparser_get_next_double(pp, -BIG_COORD, +BIG_COORD);
        ip[i].c.c[2] = argparser_get_next_double(pp, -BIG_COORD, +BIG_COORD);
        fprintf(stderr, "  ip[%d] = [ %10.2f %10.2f %10.2f ]\n", i, ip[i].c.c[0], ip[i].c.c[1], ip[i].c.c[2]);
      }
    fprintf(stderr, "\n");
    hr2_pmap_t IM = hr2_pmap_from_four_points(&(ip[0]), &(ip[1]), &(ip[2]), &(ip[3]));
    /* Parse the four output points {op[0..3]}: */
    hr2_point_t op[4];
    for (i = 0; i < 4; i++)
      { op[i].c.c[0] = 1.0;
        op[i].c.c[1] = argparser_get_next_double(pp, -BIG_COORD, +BIG_COORD);
        op[i].c.c[2] = argparser_get_next_double(pp, -BIG_COORD, +BIG_COORD);
        fprintf(stderr, "  op[%d] = [ %10.2f %10.2f %10.2f ]\n", i, op[i].c.c[0], op[i].c.c[1], op[i].c.c[2]);
      }
    fprintf(stderr, "\n");
    hr2_pmap_t OM = hr2_pmap_from_four_points(&(op[0]), &(op[1]), &(op[2]), &(op[3]));
    /* Compute the projective map: */
    hr2_pmap_t M;
    r3x3_mul(&(IM.inv), &(OM.dir), &(M.dir));
    r3x3_mul(&(OM.inv), &(IM.dir), &(M.inv));
    return M;
  }
    
hr2_pmap_t argparser_get_proj_map(argparser_t *pp)
  {
    hr2_pmap_t M;
    if (argparser_keyword_present(pp, "-matrix"))
      { M = argparser_get_next_proj_map_matrix(pp); }
    else if (argparser_keyword_present(pp, "-points"))
      { M = argparser_get_next_proj_map_from_points(pp); }
    else 
      { r3x3_ident(&(M.dir));  r3x3_ident(&(M.inv)); }
    return M;
  }
    
hr2_pmap_t argparser_get_next_feature_map(argparser_t *pp)
  {
    r3x3_t M;
    r3x3_ident(&M);

    for (int32_t j = 1; j < 3; j++)
      { M.c[0][j] = argparser_get_next_double(pp, -1000.0, +1000.0); }
    
    for (int32_t i = 1; i < 3; i++)
      { for (int32_t j = 1; j < 3; j++)
          { M.c[i][j] = argparser_get_next_double(pp, -1000.0, +1000.0); }
      }
      
    double mag = argparser_get_next_double(pp, 0.001, 1000.0);
    
    for (int32_t i = 1; i < 3; i++)
      { for (int32_t j = 1; j < 3; j++)
          { M.c[i][j] /= mag; }
      }

    hr2_pmap_t A;
    A.dir = M;
    r3x3_inv(&M, &(A.inv));
    return A;
  }

/* Copyright © 2003 by Jorge Stolfi.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appears in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty of any kind.
*/
 

/* See {obj_file_write.h}. */
/* Last edited on 2024-12-05 10:39:21 by stolfi */
 
#define obj_file_write_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <affirm.h>

#include <obj_file.h>
#include <obj_file_write.h>

#define debug FALSE

void obj_file_write_coords_table(FILE *wr, char *elname, char *cmd, r3_vec_t *P, string_vec_t *VL, int32_t prc);
  /* Writes to {wr} the 3-vectors {P.e[0..P.ne-1]}, one per line,
    preceded by the {cmd} string. Each coordinate is written with {prc}
    decimal fraction digits.  If {VL} is not null and {VL.e[k]} is not {NULL},
    also writes the string {VL.e[k]} as a '#'-comment after the 
    coordinates of each {P.e[k]}.
    
    The {elname} is used in a comment line before the data. */

void obj_file_write_face_tables
  ( FILE *wr,
    int32_t nv,
    obj_file_face_vec_t *FV,
    int32_t nt,
    obj_file_face_vec_t *FT,
    int32_t nn,
    obj_file_face_vec_t *FN
  );
  /* Writes tp {wr} the face corner info from the tables {FV,FT,FN}.
    Assumes that the definitions of {nv} vertices,{nt} textpoints, and {nn} normals have been 
    written.
    
    Each index in the tables is zero-based, whereas OBJ indices are 1-based.
    Therefore, each index {ix} from those tables is written as {ix+1}. */ 

void obj_file_write(FILE *wr, obj_file_data_t *D, int prec)
  {
    int32_t nv = D->V.ne;  /* Number of vertices */
    int32_t nt = D->T.ne;  /* Number of texpoints */
    int32_t nn = D->N.ne;  /* Number of normals */
    assert(nv == D->VL.ne);
    
    fprintf(wr, "# written by %s from JSLIBS/libquad\n", __FUNCTION__);
    fprintf(wr, "\n");
    obj_file_write_coords_table(wr, "vertex", "v", &(D->V), &(D->VL), prec);
    obj_file_write_coords_table(wr, "texpoint", "vt", &(D->T), NULL, prec);
    obj_file_write_coords_table(wr, "normal", "vn", &(D->N), NULL, 7);
    fprintf(wr, "# faces\n");
    obj_file_write_face_tables(wr, nv, &(D->FV), nt, &(D->FT), nn, &(D->FN));
    fflush(wr);
    return;
  }
    
void obj_file_write_coords_table(FILE *wr, char *elname, char *cmd, r3_vec_t *P, string_vec_t *PL, int32_t prc)
  { int32_t np = P->ne;  /* Number of entries */
    fprintf(wr, "# %s coordinates\n", elname);
    for (uint32_t kp = 0;  kp < np; kp++)
      { fprintf(wr, "%s", cmd);
        r3_t *Pk = &(P->e[kp]);
        for (uint32_t j = 0;  j < 3; j++)
          { fprintf(wr, " %.*f", prc, Pk->c[j]); }
        if (PL != NULL)
          { char *PLk = PL->e[kp];
            if (PLk != NULL) { fprintf(wr, " # %s\n", PLk); }
          }
        fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
  }

void obj_file_write_face_tables
  ( FILE *wr,
    int32_t nv,
    obj_file_face_vec_t *FV,
    int32_t nt,
    obj_file_face_vec_t *FT,
    int32_t nn,
    obj_file_face_vec_t *FN
  )
  { int32_t nf = FV->ne; /* Number of faces. */
    demand(FT->ne == nf, "inconsistent face count (FT)");
    demand(FN->ne == nf, "inconsistent face count (FN)");
    fprintf(wr, "# face corners\n");
    for (uint32_t kf = 0;  kf < nf; kf++)
      { if (debug) { fprintf(stderr, "  writing face %d\n", kf); }
        fprintf(wr, "f");
        int32_vec_t *FVk = &(FV->e[kf]);
        int32_vec_t *FTk = &(FT->e[kf]);
        int32_vec_t *FNk = &(FN->e[kf]);
        int32_t nc = FVk->ne;
        demand(nc == FTk->ne, "inconsistent corner count (FT.e[kf])");
        demand(nc == FNk->ne, "inconsistent corner count (FN.e[kf])");
        for (uint32_t kc = 0;  kc < nc; kc++)
          { int32_t ixv = FVk->e[kc]; if (ixv != -1) { ixv += 1; }
            int32_t ixt = FTk->e[kc]; if (ixt != -1) { ixt += 1; }
            int32_t ixn = FNk->e[kc]; if (ixn != -1) { ixn += 1; }
            if (debug) { fprintf(stderr, "    writing corner %d = %d/%d/%d\n", kc, ixv, ixt, ixn); }
            demand((ixv >= 1) && (ixv <= nv), "bad vertex index");
            fprintf(wr, " %d", ixv);
            if ((ixt >= 0) || (ixn >= 0))
              { fprintf(wr, "/");
                if (ixt >= 0)
                  { demand((ixt >= 1) && (ixt <= nt), "bad texpoint index");
                    fprintf(wr, "%d", ixt);
                  }
                if (ixn >= 0)
                  { demand((ixn >= 1) && (ixn <= nn), "bad normal index");
                    fprintf(wr, "/");
                    fprintf(wr, "%d", ixn);
                  }
              }
          } 
        fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
  }

/* See {obj_file.h}. */
/* Last edited on 2025-01-09 23:53:47 by stolfi */
 
#define obj_file_C_copyright \
  "Copyright © 2024 State University of Campinas (UNICAMP).\n\n" jslibs_copyright
 
/* Written by J. Stolfi in June 2024. */ 

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <affirm.h>
#include <vec.h>

#include <obj_file.h>

obj_file_data_t *obj_file_data_new(void)
  {
    obj_file_data_t *D = talloc(1, obj_file_data_t);
    D->V = r3_vec_new(100);
    D->T = r3_vec_new(100);
    D->N = r3_vec_new(100);
    D->VL = string_vec_new(100);
    uint32_t NF_alloc = 50;
    D->FV = obj_file_face_vec_new(NF_alloc);
    D->FT = obj_file_face_vec_new(NF_alloc);
    D->FN = obj_file_face_vec_new(NF_alloc);
    return D;
  }
  
bool_t obj_file_data_compare(obj_file_data_t *D1, obj_file_data_t *D2, double tol, bool_t verbose)
  { bool_t eq = TRUE; /* Until proven guilty. */
    
    auto void numcmp(char *name, int32_t n1, int32_t n2, int32_t kf, int32_t kc);
      /* Compares integers {n1} and {n2}. If they differ, sets {res}
        to {FALSE} and prints an error message assuming those
        numbers are called {name}. If {kf} is non-negative
        assumes they refer to face {kf} If {kc} is non-negative,
        assumes they refer to corner {kc} of that face. */
    
    auto void ptscmp(char *name, r3_vec_t *P1, r3_vec_t *P2, double tol);
      /* Compares point vectors {P1} and {P2}. If they differ on any coordinate 
        by more than{tol}, sets {eq} to {FALSE} and prints an error message assuming those
        vectors are called {name}. */
    
    auto void labcmp(char *name, string_vec_t *L1, string_vec_t *L2);
      /* Compares string vectors {L1} and {L2}. If they differ,
        sets {eq} to {FALSE} and prints an error message assuming those
        vectors are called {name}. */
    
    numcmp("vertex counts",   (int32_t)D1->V.ne, (int32_t)D2->V.ne, -1, -1);
    numcmp("texpoint counts", (int32_t)D1->T.ne, (int32_t)D2->T.ne, -1, -1);
    numcmp("normal counts",   (int32_t)D1->N.ne, (int32_t)D2->N.ne, -1, -1);
    assert(D1->VL.ne == D1->V.ne);

    ptscmp("V", &(D1->V), &(D2->V), tol);
    ptscmp("T", &(D1->T), &(D2->T), tol);
    ptscmp("N", &(D1->N), &(D2->N), obj_file_data_tol_normal);
    labcmp("VL", &(D1->VL), &(D2->VL));
    
    numcmp("face count",     (int32_t)D1->FV.ne, (int32_t)D2->FV.ne, -1, -1);
    assert(D1->FT.ne == D1->FV.ne); assert(D2->FT.ne == D2->FV.ne);
    assert(D1->FN.ne == D1->FV.ne); assert(D2->FN.ne == D2->FV.ne);
    if (D1->FV.ne == D2->FV.ne)
      { uint32_t NF = D1->FV.ne; /* Number of faces. */
        for (int32_t kf = 0;  kf < NF; kf++)
          { int32_vec_t *FV1k = &(D1->FV.e[kf]); int32_vec_t *FV2k = &(D2->FV.e[kf]);
            int32_vec_t *FT1k = &(D1->FT.e[kf]); int32_vec_t *FT2k = &(D2->FT.e[kf]);
            int32_vec_t *FN1k = &(D1->FN.e[kf]); int32_vec_t *FN2k = &(D2->FN.e[kf]);
            numcmp("corner counts", (int32_t)FV1k->ne, (int32_t)FV2k->ne, kf, -1);
            assert(FT1k->ne == FV1k->ne); assert(FT2k->ne == FV2k->ne); 
            assert(FN1k->ne == FV1k->ne); assert(FN2k->ne == FV2k->ne); 
            if (FV1k->ne == FV2k->ne)
              { uint32_t NC = FV1k->ne; /* Number of corners. */
                for (int32_t kc = 0;  kc < NC; kc++)
                  { numcmp("vertex indices",   FV1k->e[kc], FV2k->e[kc], kf, kc);
                    numcmp("texpoint indices", FT1k->e[kc], FT2k->e[kc], kf, kc);
                    numcmp("normal indices",   FN1k->e[kc], FN2k->e[kc], kf, kc);
                  }
              }
          }
      }
      
   return eq;
      
   void numcmp(char *name, int32_t n1, int32_t n2, int32_t kf, int32_t kc)
     { if (n1 != n2)
         { eq = FALSE;
           if (verbose) 
             { fprintf(stderr, "  !! %s", name);
               if (kc >= 0) { fprintf(stderr, " of corner %d", kc); }
               if (kf >= 0) { fprintf(stderr, " of face %d", kf); }
               fprintf(stderr, " differ: %d %d\n", n1, n2);
             }
         }
     }
   
   void ptscmp(char *name, r3_vec_t *P1, r3_vec_t *P2, double tol)
     { assert(P1->ne == P2->ne);
       uint32_t np = P1->ne;
       for (int32_t kp = 0;  kp < np; kp++)
         { r3_t *P1k = &(P1->e[kp]);
           r3_t *P2k = &(P2->e[kp]);
           for (int32_t ka = 0;  ka < 3; ka++)
             { double P1c = P1k->c[ka];
               double P2c = P2k->c[ka];
               if (fabs(P1c - P2c) >tol)
                 { eq = FALSE;
                   if (verbose)
                     { fprintf(stderr, "  !! coordinates %d of %s[%d] differ:", ka, name, kp);
                       fprintf(stderr, " %.7f %.7f\n", P1c, P2c);
                     }
                 }
             }
         }
     }

   void labcmp(char *name, string_vec_t *L1, string_vec_t *L2)
     { assert(L1->ne == L2->ne);
       uint32_t NS = L1->ne;
       for (int32_t ks = 0;  ks < NS; ks++)
         { char *L1k = L1->e[ks];
           char *L2k = L2->e[ks];
           if (strcmp(L1k, L2k) != 0)
             { eq = FALSE;
               if (verbose)
                 { fprintf(stderr, "  !! labels of %s[%d] differ:", name, ks);
                   fprintf(stderr, " \"%s\" \"%s\"\n", L1k, L2k);
                 }
             }
         }
     }

  }

void obj_file_data_free(obj_file_data_t *D)
  {
    free(D->V.e);
    free(D->T.e);
    free(D->N.e);
    uint32_t NF = D->FV.ne;
    assert(D->FT.ne == NF);
    assert(D->FN.ne == NF);
    for (int32_t kf = 0;  kf < NF; kf++)
      { free(D->FV.e[kf].e); 
        free(D->FT.e[kf].e); 
        free(D->FN.e[kf].e); 
      }
    free(D->FV.e);
    free(D->FT.e);
    free(D->FN.e);
    free(D);
  }
 
vec_typeimpl(obj_file_face_vec_t,obj_file_face_vec,int32_vec_t);


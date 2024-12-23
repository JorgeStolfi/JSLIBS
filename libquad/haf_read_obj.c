/* See {haf_read_obj.h}. */
/* Last edited on 2024-12-22 10:28:29 by stolfi */
 
#define haf_read_obj_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright
 
/* Written by J. Stolfi in June 2024. */ 

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <r3.h>
#include <jsstring.h>
#include <jsprintf.h>
#include <affirm.h>
#include <fget.h>

#include <obj_file.h>
#include <obj_file_read.h>

#include <haf.h>
#include <haf_read_obj.h>

#define debug FALSE

void haf_read_obj_err(const char *func, char *msg);
  /* Prints "{func}: ** error {msg}\n" to {stderr} and bombs out. */

void haf_read_obj_file
  ( FILE *rd,
    haf_edge_id_t eid0,
    uint32_t *nf_P,
    r3_vec_t *V_P,
    string_vec_t *VL_P,
    haf_arc_vec_t *A_P,
    int32_t **fleft_P,
    int32_t **vorg_P,
    bool_t verbose
  )
  { 
    const char *func = __FUNCTION__; /* For error messages. */ 
    
    obj_file_data_t *D = obj_file_read(rd, verbose);
    
    uint32_t nf = D->FV.ne;          /* Number of faces ('f' lines) found. */
    uint32_t nv = D->V.ne;       /* Number of vertices ('v' lines) specified in the OBJ file. */
    assert(D->VL.ne == nv);
    
    /* Determine the total number {ns_max} of sides in the faces: */
    uint32_t ns_max = 0;
    for (int32_t kf = 0;  kf < nf; kf++)
      { int32_vec_t *FVk = &(D->FV.e[kf]);
        uint32_t nck = FVk->ne;
        demand(nck >= 3, "face has fewer than 3 corners");
        ns_max += nck;
      }
    if (verbose) { fprintf(stderr, "input mesh has %d total face sides\n", ns_max); }
    
    /* Alocate the edge table: */
    uint32_t ne_max = ns_max; /* In case there is free border. */
    if (verbose) { fprintf(stderr, "allocating the edge table {A} for up to %d edges...\n", ne_max); }
    haf_arc_vec_t A = haf_arc_vec_new(ne_max);
    for (int32_t ke = 0;  ke < ne_max; ke++) { A.e[ke] = NULL;  }
      
    /* Alocate the arc tables: */
    uint32_t na_max = 2*ne_max; 
    if (verbose) { fprintf(stderr, "allocating the {vorg,fleft} tables table for up to %d arcs...\n", na_max); }
    int32_t *fleft = talloc(na_max, int32_t);  /* In case there is free border. */
    int32_t *vorg = talloc(na_max, int32_t);   /* In case there is free border. */

    for (int32_t ka = 0;  ka < na_max; ka++) { fleft[ka] = -1; vorg[ka] = -1; }
        
    /* Vertex tables: */
    /* {V,VL} are copies of {D->V,D->VL}: */
    if (verbose) { fprintf(stderr, "copying the {V} and {VL} tables (%d entries)...\n", nv); }
    r3_vec_t V = r3_vec_new(nv);
    string_vec_t VL = string_vec_new(nv);
    /* {vdeg[kv]} is out-degree of vertex {D->V.e[kv]}: */
    uint32_t *vdeg = talloc(nv, uint32_t);  
    /* {aout[kv].e[0..vdeg[kv]-1]} are the arcs out of {D->V.e[kv]}: */
    haf_arc_vec_t *aout = talloc(nv, haf_arc_vec_t);

    for (int32_t kv = 0;  kv < nv; kv++)
      { V.e[kv] = D->V.e[kv];
        /* Make a copy of the label: */
        char *lab = D->VL.e[kv];
        VL.e[kv] = (lab == NULL ? NULL : txtcat(lab,"")); 
        vdeg[kv] = 0;
        aout[kv] = haf_arc_vec_new(5);  /* To be expanded as needed. */
      }
    
    uint32_t ne = 0; /* Number of edge records actually created. */

    auto void add_side(int32_t kv_org, int32_t kv_dst, int32_t kf);
      /* Checks {aout[kv_org}} and {aout[kv_dst]} to see whether
        the unordered pair {{kv_org,kv_dst}} is a new edge. If so, 
        creates its record with ID {ne}, adds
        its base arc (which will have origin {kv_org}) {A} as {A[ne]} and
        increments {ne}. Also inserts both arcs in
        {aout}, and sets {vorg[ka],vorg[kb]} appropriately.  Otherwise,
        if the edge already exists, checks that its base arc has origin {kv_dst},
        and is stored in {A} at the proper index. In any case
        sets {fleft[ka]} to {kf}. */
    
    auto void get_arc_indices
      ( haf_arc_t r,
        int32_t *kr_P,
        haf_arc_t s,
        int32_t *ks_P,
        int32_t *ke_P,
        haf_edge_id_t ide_min,
        haf_edge_id_t ide_max
      );
      /* Obtains and checks the indices {kr,ks} of arcs {r} and {s},
        and the index {ke} of their edge, from their indices.
        Arc {r} must be the base arc of the edge. Also 
        checks that the edge ID is in the range {ide_min..ide_max}. */
    
    auto haf_arc_t find_in_aout(int32_t kv0, int32_t kv1, int32_t kf);
      /* Checks {aout[kv0]} to see whether the arc {(kv0-> kv1)} 
        has already been used. If so, returns that arc.
        Otherwise, returns {NULL}.  */
      
    auto void check_not_in_aout(int32_t kv0, int32_t kv1, int32_t kf);
      /* Checks whether the arc {(kv0-> kv1)} 
        has already been used.  If so, aborts with an error message. */

    auto void prtvert(char *pref, int32_t kv, char *suff);
      /* Prints the index {kv} ov vertex {V.e[kv]}, and its label {VL.e[kv]}. */
      
    auto void prtarc(char *pref, haf_arc_t a, char *suff);
      /* Prints the data of arc {a}, including the edge index and ID,
        the direction bit, and the origin and destination vertices. */

    if (verbose) { fprintf(stderr, "scanning the faces and creating the edge records...\n"); }
    for (int32_t kf = 0;  kf < nf; kf++)
      { int32_vec_t *FVk = &(D->FV.e[kf]);
        uint32_t nck = FVk->ne; /* Number of arcs in face border. */
        if (debug) { fprintf(stderr, "  processing face F[%d] (%d corners)\n", kf, nck); }
        /* Scan the border to create the edge records and set arc tables: */
        int32_t kv_prev = FVk->e[nck-1];  assert((kv_prev >= 0) && (kv_prev < D->V.ne));
        for (int32_t kc = 0;  kc < nck; kc++)
          { int32_t kv_this = FVk->e[kc]; assert((kv_this >= 0) && (kv_this < D->V.ne));
            add_side((int32_t)kv_prev, (int32_t)kv_this, kf);
            kv_prev = kv_this;
          }
        /* Scan the border again to set the {next} links: */
        int32_t kv0_ini = FVk->e[nck-1]; assert((kv0_ini >= 0) && (kv0_ini < D->V.ne));
        int32_t kv1_ini = FVk->e[0]; assert((kv1_ini >= 0) && (kv1_ini < D->V.ne));
        haf_arc_t a_prev = find_in_aout((int32_t)kv0_ini, (int32_t)kv1_ini, kf);
        assert(a_prev != NULL);
        for (int32_t kc = 0;  kc < nck; kc++)
          { int32_t kv0 = FVk->e[kc]; assert((kv0 >= 0) && (kv0 < D->V.ne));
            int32_t kv1 = FVk->e[(kc+1)%(int32_t)nck]; assert((kv1 >= 0) && (kv1 < D->V.ne));
            haf_arc_t a_this = find_in_aout((int32_t)kv0, (int32_t)kv1, kf);
            assert(a_this != NULL);
            haf_set_lnext(a_prev, a_this);
            a_prev = a_this;
          }
      }
    haf_arc_vec_trim(&A, ne);
    if (verbose) { fprintf(stderr, "found %d edges\n", ne); }
   
    uint32_t na = 2*ne;
    if (na != na_max)
      { assert(na < na_max);
        vorg = retalloc(vorg, na, int32_t);
        fleft = retalloc(fleft, na, int32_t);
      }
     
    /* Check for free borders: */
    if (verbose) { fprintf(stderr, "checking for free borders...\n"); }
    uint32_t na_free = 0;  /* Number of unused arc (free border edges). */
    for (int32_t ka = 0;  ka < 2*ne; ka++)
      { if (fleft[ka] == -1)
          { haf_arc_t a = A.e[ka];
            assert(a != NULL);
            haf_arc_t b = haf_sym(a);
            assert(b != NULL);
            int32_t kb = (int32_t)(haf_arc_id(a) - 2*eid0);
            int32_t kv_org = vorg[ka]; assert((kv_org >= 0) && (kv_org < D->V.ne));
            int32_t kv_dst = vorg[kb]; assert((kv_dst >= 0) && (kv_dst < D->V.ne));
            fprintf(stderr, "  arc %d = (%d->%d) has only one face\n", ka, kv_org,kv_dst);
            na_free++;
          }
      }
    if (na_free > 0) { fprintf(stderr, "mesh has %d free border edges\n", na_free); }
         
    /* We no longer need {aout}: */
    if (verbose) { fprintf(stderr, "discarding the {aout} table...\n"); }
    for (int32_t kv = 0;  kv < nv; kv++) { free(aout[kv].e); }
    free(aout); aout = NULL;
    
    /* Check for unused vertices and squeeze them out: */
    if (verbose) { fprintf(stderr, "deleting unused vertices...\n"); }
    int32_t *kv_map = talloc(nv, int32_t); /* Maps original vert indices to squeezed ones. */
    uint32_t nv_old = nv; /* Number of vertices in OBJ file. */
    uint32_t nv_new = 0;  /* Number of vertices retained. */
    for (int32_t kv = 0;  kv < nv; kv++)
      { if (vdeg[kv] == 0)
          { if (verbose) { fprintf(stderr, "  vertex V[%d] of file is not used, discarded\n", kv); }
            free(VL.e[kv]); 
            VL.e[kv] = NULL;
            kv_map[kv] = -1;
          }
        else
          { /* Vertex is used: */
            if (debug) { fprintf(stderr, "  vertex %d of file is is renumbered to %d\n", kv, nv_new); }
            if (nv_new < kv)
              { V.e[nv_new] = V.e[kv];
                VL.e[nv_new] = VL.e[kv]; 
                VL.e[kv] = NULL;
              }
            kv_map[kv] = (int32_t)nv_new; 
            nv_new++;
          }
      }
    r3_vec_trim(&V, nv_new);
    string_vec_trim(&VL, nv_new);
    if (verbose) { fprintf(stderr, "found %d used vertices out of %d\n", nv_new, nv); }
    
    /* Map {vorg} to account for vertex squeezing: */
    if (nv_new < nv)
      { if (verbose) { fprintf(stderr, "renumbering used vertices\n"); }
        for (int32_t ka = 0;  ka < na; ka++)
          { int32_t kv_old = vorg[ka];
            if (debug) { fprintf(stderr, "    vorg[%d] = V[%d] fleft[%d] = F[%d]", ka, kv_old, ka, fleft[ka]); } 
            assert((kv_old >= 0) && (kv_old < nv_old));
            int32_t kv_new = kv_map[kv_old];
            assert((kv_new >= 0) && (kv_new < nv_new));
            if (debug) { fprintf(stderr, " -> V[%d](%s)\n", kv_new, VL.e[kv_new]); } 
            vorg[ka] = kv_new;
          }
      }
    nv = nv_new;
  
    /* We don't need these any more: */
    if (verbose) { fprintf(stderr, "discarding work tables...\n"); }
    free(vdeg);
    free(kv_map);
    
    if (verbose) { fprintf(stderr, "discarding the {obj_file_data_t} object...\n"); }
    obj_file_data_free(D);
    
    /* Return results: */
    (*nf_P) = nf;
    (*A_P) = A;
    (*V_P) = V;
    (*VL_P) = VL;
    (*fleft_P) = fleft;
    (*vorg_P) = vorg;
    
    return;
    
    void add_side(int32_t kv_org, int32_t kv_dst, int32_t kf)
      { 
        if (debug) { fprintf(stderr, "    adding arc V[%d] -> V[%d] on border of F[%d]\n", kv_org, kv_dst, kf); }
        if ((kv_org < 0) || (kv_org >= nv))
          { 
            char *msg = jsprintf("corner of F[%d] is an invalid vertex V[%d]", kf, kv_org);
            haf_read_obj_err(func, msg);
          }
        if ((kv_dst < 0) || (kv_dst >= nv))
          { 
            char *msg = jsprintf("corner of F[%d] is an invalid vertex V[%d]", kf, kv_dst);
            haf_read_obj_err(func, msg);
          }
        if (kv_org == kv_dst)
          { 
            char *msg = jsprintf("arc on face F[%d] is a loop V[%d]--V[%d]", kf, kv_org, kv_dst);
            haf_read_obj_err(func, msg);
          }
        /* Checks {aout[kv_org]} if  this arc {kv_out->kv_dst} occurred before: */
        check_not_in_aout(kv_org, kv_dst, kf);
        haf_arc_t a = NULL;
        haf_arc_t b = find_in_aout(kv_dst, kv_org, kf);
        int32_t ka = -1, kb = -1, ke = -1;
        if (b == NULL)
          { /* Edge is new, create: */
            if (debug) { fprintf(stderr, "    new edge A[%d]\n", ne); }
            a = haf_make_stick(ne + eid0);
            b = haf_sym(a);
            /* Get the indices of {a} and {b} (note {a} is the base arc): */
            ke = (int32_t)ne;
            ka = (int32_t)(2*ne);
            kb = ka + 1;
            /* Save {a} in {A} table: */
            assert(A.e[ke] == NULL);
            A.e[ke] = a;
            /* Set the {vorg} indices of {a,b}: */
            assert(vorg[ka] == -1);
            vorg[ka] = kv_org;
            assert(vorg[kb] == -1);
            vorg[kb] = kv_dst;
            /* We'll set {fleft} of {a} right down below, but we don't know the {fleft} of {b} yet: */
            assert(fleft[kb] == -1);
            ne++;
          }
        else
          { /* Edge was created previously when its {haf_sym} was found: */
            a = haf_sym(b);
            /* Get indices (IDs) of {a} and {b} (note {b} is the base arc): */
            get_arc_indices(b, &kb, a, &ka, &ke, eid0, eid0+ne-1);
            if (debug) { prtarc("    twin of old edge ", b,  ""); fprintf(stderr, " from face F[%d]\n", fleft[kb]); }
            /* Most table entries for {a} must have been set: */
            assert(A.e[ke] == b);
            assert(vorg[ka] == kv_org);
            assert(vorg[kb] == kv_dst);
            assert(fleft[kb] != -1);
          }
        /* Set the {fleft} index of {a}:*/
        assert(fleft[ka] == -1);
        fleft[ka] = kf;
        
        /* Save arc {a} in {aout[kv_org]}: */ 
        haf_arc_vec_t *aout_org = &(aout[kv_org]);
        int32_t ja = (int32_t)vdeg[kv_org];
        haf_arc_vec_expand(aout_org, ja);
        aout_org->e[ja] = a;
        vdeg[kv_org]++;
        assert(vorg[ka] == kv_org);
        assert(vorg[kb] == kv_dst);
        if (debug) { prtarc("    added arc ", a, "\n"); }
      }
          
    void get_arc_indices
      ( haf_arc_t r,
        int32_t *kr_P,
        haf_arc_t s,
        int32_t *ks_P,
        int32_t *ke_P,
        haf_edge_id_t ide_min,
        haf_edge_id_t ide_max
      )
      { 
        if (debug) 
          { fprintf(stderr, "      getting arc indices"),
            fprintf(stderr, " r = %lu = %lu:%u", haf_arc_id(r), haf_edge_id(r), haf_dir_bit(r));
            fprintf(stderr, " s = %lu = %lu:%u", haf_arc_id(s), haf_edge_id(s), haf_dir_bit(s));
            fprintf(stderr, " e = %lu = A[%d]", haf_edge_id(r), (int32_t)(haf_edge_id(r) - eid0) );
            fprintf(stderr, "\n");
          }
        /* Get indices (IDs) of {r} and {s}: */
        haf_edge_id_t ide = haf_edge_id(r);
        assert((ide >= ide_min) && (ide <= ide_max));
        /* Paranoia: */
        haf_arc_id_t ida = haf_arc_id(r); 
        assert(ida == 2*ide);
        haf_arc_id_t idb = haf_arc_id(s);
        assert(idb == ida + 1);
        /* Convert {ID}s to indices: */
        int32_t ke = (int32_t)(ide - eid0);
        int32_t kr = 2*ke;
        int32_t ks = 2*ke+1;
        /* Return: */
        (*kr_P) = kr; (*ks_P) = ks; (*ke_P) = ke;
      }

    haf_arc_t find_in_aout(int32_t kv0, int32_t kv1, int32_t kf)
      { /* Checks {aout[kv0]} if arc {kv0->kv1} occurred before: */
        haf_arc_vec_t *aout0 = &(aout[kv0]);
        for (int32_t ja = 0;  ja < vdeg[kv0]; ja++)
          { haf_arc_t a = aout0->e[ja];
            assert (a != NULL);  /* We should never store a {NULL} there. */
            int32_t ka = (int32_t)(haf_arc_id(a) - 2*eid0);
            int32_t kv_org_a = vorg[ka];
            assert(kv_org_a == kv0); /* We should only store arcs out of {kv0} there. */
            haf_arc_t b = haf_sym(a);
            assert (b != NULL);  /* Axiom. */
            int32_t kb = (int32_t)(haf_arc_id(b) - 2*eid0);
            int32_t kv_dst_a = vorg[kb];
            if (kv_dst_a == kv1) { return a; }
          }
        return NULL;
     }
    
    void check_not_in_aout(int32_t kv0, int32_t kv1, int32_t kf)
      { haf_arc_t a = find_in_aout(kv0, kv1, kf);
        if (a != NULL)
          { /* This arc has been used before: */
            int32_t ka = (int32_t)(haf_arc_id(a) - 2*eid0);
            int32_t kb = (int32_t)(haf_arc_id(haf_sym(a)) - 2*eid0);
            int32_t kf_old = fleft[ka];
            assert(kf_old != -1);
            assert(vorg[ka] == kv0);
            assert(vorg[kb] == kv1);
            prtarc("arc ", a, "\n");
            char *msg = jsprintf("arc of face F[%d] above also used in face F[%d]", kf, kf_old);
            haf_read_obj_err(func, msg);
          }
        return;
      }
          
    void prtarc(char *pref, haf_arc_t a, char *suff)
      { haf_arc_t b = haf_sym(a);
        int32_t ka = (int32_t)(haf_arc_id(a) - 2*eid0);
        int32_t kb = (int32_t)(haf_arc_id(b) - 2*eid0);
        haf_edge_id_t eid = haf_edge_id(a);
        haf_dir_bit_t adir = haf_dir_bit(a);
        int32_t ke = (int32_t)(eid - eid0);
        char *op = (adir == 0 ? "" : ".sym");
        assert(A.e[ke] == haf_base_arc(a));
        fprintf(stderr, "%sA[%d]%s(%lu:%u)", pref, ke, op, eid, adir);
        prtvert(" = ", vorg[ka], "");
        prtvert(" --> ", vorg[kb], "");
        int32_t kf = fleft[ka];
        if (kf != -1) { fprintf(stderr, " of face F[%d]", kf); }
        fprintf(stderr, "%s", suff);
      }
      
    void prtvert(char *pref, int32_t kv, char *suff)
      { fprintf(stderr, "%s", pref);
        fprintf(stderr, "V[%d]", kv);
        if ((kv >= 0) && (kv < nv))
          { char *lab = VL.e[kv];
            fprintf(stderr, "(%s)", (lab == NULL ? "" : lab));
          }
        fprintf(stderr, "%s", suff);
      }
  }

void haf_read_obj_err(const char *func, char *msg)
  { fprintf(stderr, "%s: ** %s\n", func, msg);
    assert(FALSE);
  }

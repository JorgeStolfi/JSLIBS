/* See {tosl_build_lists_hash.h} */
/* Last edited on 2024-10-06 16:50:14 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_arc_list.h>
#include <tosl_build_lists_hash.h>

int32_t tosl_build_lists_hash_func(tosl_coord_t  Z, int32_t NH, tosl_coord_t  Zmin, tosl_coord_t  Zmax);
  /* Hashes a {Z}-cordinate in {Zmin..Zmax} to an index in {0..NH-1}. A monotonic non-decreasing function of {Z}. */

tosl_plane_id_t *tosl_build_lists_hash_make_table(int32_t NH, int32_t NP, tosl_coord_t Zplane[], int32_t debug);
  /* Builds a hash table {iphahs[0..NH-1]} such that, for any coordinate {Z},
    {iphash[hash(Z)]} is an index {ip} such that {Zpane[ip]} is close enough 
    to {Z}. */

tosl_arc_id_t *tosl_build_lists_hash
  ( int32_t NP,
    tosl_coord_t Zplane[],
    tosl_mesh_t *mesh,
    int32_t debug
  )
  {
    assert((mesh->NE > 0) && (mesh->NE <= tosl_mesh_MAX_EDGES));
    assert((NP > 0) && (NP <= tosl_mesh_MAX_PLANES));

    int32_t NA = 2*mesh->NE;
    
    /* Allocate the bucket table and set all entries to empty:  */
    tosl_arc_id_t *L = malloc(NP*sizeof(tosl_arc_id_t));
    assert(L != NULL);
    for (size_t ip = 0; ip < NP; ip++) { L[ip] = -1; }
    
    /* Hash table for {Zplane[0..NP-1]}: */
    int32_t NH = 2*NP;
    tosl_plane_id_t *iphash = tosl_build_lists_hash_make_table(NH, NP, Zplane, debug);
    
    int32_t ns = 0; /* Number of sequential search steps performed. */
    for (tosl_arc_id_t ia = 0; ia < NA; ia++)
      { /* Make {Arc[ia]} into a singleton list: */
        mesh->Arc[ia].pred = ia;
        mesh->Arc[ia].succ = ia;
        
        /* Get the origin and destination {Z}-coordinates {Zorg,Zdst}: */
        tosl_coord_t Zorg = mesh->Vpos[mesh->Arc[ia].ivorg].c[2]; /* {Z} of origin of {ia}.  */
        tosl_arc_id_t ja = tosl_sym(ia);
        tosl_coord_t Zdst = mesh->Vpos[mesh->Arc[ja].ivorg].c[2]; /* {Z} of destination of {ia}.  */
        
        /* Add the arc to a list {L[ip]} if appropriate: */
        if ((Zorg < Zdst) && (Zorg < Zplane[NP-1]) && (Zdst > Zplane[0]))
          { /* Arc is up-going and is not above or below all planes. */
            tosl_plane_id_t ip = -1;
            if (Zorg < Zplane[0]) 
              { ip = 0; }
            else
              { int32_t ih = tosl_build_lists_hash_func(Zorg, NH, Zplane[0], Zplane[NP-1]);
                assert ((ih >= 0) && (ih < NH));
                ip = iphash[ih];
                if (debug) { fprintf(stderr, "  hash: %8d", ip); }
                assert((ip >= 0) && (ip < NP));
                while ((ip < NP) && (Zorg > Zplane[ip])) { ip++; ns++; }
                if (debug) { fprintf(stderr, " -> %8d", ip); }
                assert(ip < NP); /* Because {Zorg <= Zplane[NP]} */
                while ((ip > 0) && (Zorg < Zplane[ip-1])) { ip--; ns++; }
                assert(ip > 0); /* Because {Zorg >= Zplane[0]}. */
                if (debug) { fprintf(stderr, " -> %8d\n", ip); }
                assert(Zorg != Zplane[ip-1]); /* No plane-edge collision. */
                assert(Zorg != Zplane[ip]); /* No plane-edge collision. */
              }
            /* Paranoia, checking search correctness: */
            assert((ip >= 0) && (ip < NP));
            assert(Zorg < Zplane[ip]);
            assert((ip == 0) || (Zorg > Zplane[ip-1]));
            if (Zdst > Zplane[ip]) 
              { /* Arc {ia} crosses at least one slicing plane: */
                tosl_arc_list_add(&(L[ip]), ia, mesh);
              }
            else
              { /* Check that no vertices are on planes: */
                assert(Zdst < Zplane[ip]); 
              }
           }
        else
          { /* No plane-edge collisions: */
            assert((Zorg >= Zdst) || (Zorg > Zplane[NP-1]) || (Zdst < Zplane[0]));
          }
      }
    if (debug)
      { double nspe = ((double)ns)/((double)mesh->NE);
        fprintf(stderr, "performed %.2f sequential steps per edge\n", nspe);
      }
    
    free(iphash);
    return L;
  }
      
int32_t tosl_build_lists_hash_func(tosl_coord_t  Z, int32_t NH, tosl_coord_t  Zmin, tosl_coord_t  Zmax)
  { /* Affine map from {Zmin..Zmax} to {0..NH-1}: */
    int32_t ih = (int32_t)floor(((double)NH)*(Z - Zmin)/((double)Zmax - Zmin + 0.0001) + 0.5);
    if (ih < 0)
      { return 0; }
    else if (ih >= NH)
      { return NH-1; }
    else
      { return ih; }
  }

tosl_plane_id_t *tosl_build_lists_hash_make_table(int32_t NH, int32_t NP, tosl_coord_t Zplane[], int32_t debug)
  { 
    tosl_plane_id_t *iphash = malloc(NH*sizeof(tosl_plane_id_t));
    assert(iphash != NULL);
    for (int32_t ih = 0; ih < NH; ih++) { iphash[ih] = -1; }
    
    tosl_coord_t Zmin = Zplane[0];
    tosl_coord_t Zmax = Zplane[NP-1];
    
    if (debug) 
      { double tau = ((double)Zmax - Zmin)/((double)NH);
        fprintf(stderr, "  building a hash table with %d entries (approx Z resolution %.2f)\n", NH, tau);
      }
    /* Hash the slicing plane indices into {iphash[0..NH-1]} by their {Z}: */
    for (tosl_plane_id_t ip = 0; ip < NP; ip++)
      { if (ip > 0) { assert(Zplane[ip] > Zplane[ip-1]); }
        int32_t ih = tosl_build_lists_hash_func(Zplane[ip], NH, Zmin, Zmax);
        iphash[ih] = ip;
      }
      
    /* Fill the empty entries of the hash table, and check for too many collisions: */
    int32_t max_dip = 0;
    tosl_plane_id_t iplast = 0;
    int32_t ih = 0;
    while (ih < NH)
      { if (iphash[ih] == -1)
          { /* Look forward for next defined entry: */
            int32_t jh = ih + 1;
            while ((jh < NH) && (iphash[jh] == -1)) { jh++; }
            if (debug && (jh - ih >= 5)) { fprintf(stderr, "    gap %d..%d  iplast = %d\n", ih, jh-1, iplast); }
            if (jh >= NH)
              { /* Fill whole gap with {iplast}: */
                for (int32_t kh = ih; kh < NH; kh++) { iphash[kh] = iplast; }
              }
            else
              { /* Split gap between {iplast} and {iphash[jh]}: */
                tosl_plane_id_t ipnext = iphash[jh];
                int32_t kh = jh - 1;
                while (ih <= kh)
                  { iphash[ih] = iplast; 
                    iphash[kh] = ipnext;
                    ih++; kh--;
                  }
              }
            /* Skip over the gap: */
            ih = jh;
          }
          
        int32_t dip = iphash[ih] - iplast;
        assert(dip >= 0);
        if (dip > max_dip) { max_dip = dip; }
        iplast = iphash[ih];
        ih++;
      }
    if (NP >= 2) 
      { assert(max_dip >= 1);
        if (debug) { fprintf(stderr, "    max dropped planes = %d\n", max_dip - 1); }
        if ((1 << (max_dip/2)) > NP)
          { fprintf(stderr, "!! warning: there are hash table entries with %d collisions.", max_dip);
            fprintf(stderr, "  Binary search might be faster.\n");
          }
      }
    return iphash;
  }
  

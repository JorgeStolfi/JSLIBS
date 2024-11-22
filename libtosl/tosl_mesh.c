/* See {tosl_mesh.h} */
/* Last edited on 2024-11-20 05:22:06 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include <tosl.h>
#include <tosl_mesh.h>

void tosl_mesh_arc_print(FILE *wr, char *pref, tosl_arc_id_t ka, char *suff, tosl_mesh_t *mesh)
  { if (pref != NULL) { fputs(pref, wr); }
    int32_t NA = 2*mesh->NE;
    assert((ka >= 0) && (ka < NA));
    tosl_arc_t *arc = &(mesh->Arc[ka]);
    char *xka = tosl_arc_id_to_string(ka);
    char *xskip = tosl_arc_id_to_string(arc->skip);
    char *xpred = tosl_arc_id_to_string(arc->pred);
    char *xsucc = tosl_arc_id_to_string(arc->succ);
    fprintf(wr, "%-8s = { .skip = %-8s", xka, xskip);
    tosl_vert_id_t kv = arc->ivorg;
    tosl_mesh_vert_print(wr, "  .ivorg = ", kv, NULL, mesh);
    fprintf(wr, "  .pred = %s  .succ = %s }", xpred, xsucc);
    free(xka); free(xskip); free(xpred); free(xsucc);
    if (suff != NULL) { fputs(suff, wr); }
    fflush(wr);
  }

void tosl_mesh_vert_print(FILE *wr, char *pref, tosl_vert_id_t kv, char *suff, tosl_mesh_t *mesh)
  {
    if (pref != NULL) { fputs(pref, wr); }
    if ((kv >= 0) && (kv < mesh->NV)) 
      { fprintf(wr, "v%-7d", kv);
        if (mesh->Vlab != NULL)
          { if (mesh->Vlab[kv] != NULL)
              { fprintf(wr, " = %-8s", mesh->Vlab[kv]); }
            else
              { fprintf(wr, "   %-8s", ""); }
          }
        tosl_point_t *v = &(mesh->Vpos[kv]);
        fprintf(wr, " = (");
        for (int32_t j = 0; j < 3; j++) { fprintf(wr, " %8d", v->c[j]); }
        fprintf(wr, " )");
      }
    else
      { fprintf(wr, "    v%-7s", "???");
        fprintf(wr, "   %-8s", "");
        fprintf(wr, "   ");
        for (int32_t j = 0; j < 3; j++) { fprintf(wr, " %8s", ""); }
        fprintf(wr, "  ");
      }
    if (suff != NULL) { fputs(suff, wr); }
    fflush(wr);
  }

void tosl_mesh_check(tosl_mesh_t *mesh)
  {
    int32_t NA = 2*mesh->NE;
    tosl_point_t vmin, vmax;
    tosl_mesh_coord_range_get(mesh, &vmin, &vmax);
    tosl_mesh_coord_range_print(stderr, "  vertex Z range = ", &vmin, &vmax, "\n");

    /* {nd[deg]} is count of faces with degree {deg}. */
    int32_t *nd = malloc((NA+1)*sizeof(int32_t)); 
    for (int32_t deg = 0; deg <= NA; deg++) { nd[deg] = 0; }
    
    /* {seen[ia]} is true if arc {ia} was checked. */
    int8_t *seen = malloc(NA*sizeof(int8_t));
    for (int32_t ia = 0; ia < NA; ia++) { seen[ia] = 0; }

    int32_t kf = 0; /* Num of faces seen. */
    for (tosl_arc_id_t ia = 0; ia < NA; ia++)
      { 
        if (seen[ia] == 0)
          { /* New face: */
            kf++;
            int32_t deg = 0;
            tosl_arc_id_t ka = ia;
            do 
              { 
                assert(seen[ka] == 0);

                tosl_arc_id_t ka1 = tosl_sym(ka);
                tosl_arc_id_t ja = mesh->Arc[ka].skip;

                if (ja < 0)
                  { fprintf(stderr, "{.skip} undefined\n");
                    tosl_mesh_arc_print(stderr, "  ka =       ", ka, "\n", mesh); 
                    tosl_mesh_arc_print(stderr, "  sym(ka) =  ", ka1, "\n", mesh);
                    assert(0);
                  }

                assert(ja < NA);

                if (mesh->Arc[ka1].ivorg != mesh->Arc[ja].ivorg)
                  { fprintf(stderr, "{.skip}/{.ivorg} inconsistency\n");
                    tosl_arc_id_t ja1 = tosl_sym(ja);
                    tosl_mesh_arc_print(stderr, "  ka =       ", ka, "\n", mesh); 
                    tosl_mesh_arc_print(stderr, "  sym(ka) =  ", ka1, "\n", mesh); 
                    tosl_mesh_arc_print(stderr, "  ja =       ", ja, "\n", mesh);
                    tosl_mesh_arc_print(stderr, "  sym(ja) =  ", ja1, "\n", mesh);
                    assert(0);
                  }

                seen[ka] = 1;
                ka = ja; 
                deg++; 
              } while (ka != ia);

            assert((deg >= 1) && (deg <= NA));
            nd[deg]++;
          }
      }
    fprintf(stderr, "found %d faces\n", kf);
    for (int32_t deg = 0; deg <= NA; deg++)
      { if (nd[deg] != 0)
          { fprintf(stderr, "  %6d faces of degree %d\n", nd[deg], deg); }
      }
    free(seen);
    free(nd);
    return;
  }

void tosl_mesh_print(FILE *wr, tosl_mesh_t *mesh)
  {
    fprintf(wr, "  mesh has %d vertices:\n", mesh->NV);
    for (tosl_vert_id_t iv = 0; iv < mesh->NV; iv++)
      { tosl_mesh_vert_print(wr, "    ", iv, "\n", mesh); }
    fprintf(wr, "\n");

    tosl_point_t vmin, vmax;
    tosl_mesh_coord_range_get(mesh, &vmin, &vmax);
    tosl_mesh_coord_range_print(wr, "  vertex coord ranges = ", &vmin, &vmax, "\n");
    fprintf(wr, "\n");

    fprintf(wr, "  mesh has %d edges:\n", mesh->NE);
    int32_t NA = 2*mesh->NE;
    for (tosl_arc_id_t ia = 0; ia < NA; ia++)
      { tosl_mesh_arc_print(wr, "    ", ia, "\n", mesh);
        if ((ia%2) == 1) { fprintf(wr, "\n"); }
      }
    fprintf(wr, "\n");

    fprintf(wr, "  mesh has the following faces:\n");
    uint8_t *seen = malloc(NA*sizeof(uint8_t));
    assert(seen != NULL);
    for (tosl_arc_id_t ia = 0; ia < NA; ia++) { seen[ia] = 0; }

    int32_t kf = 0; /* Faces found so far. */
    for (tosl_arc_id_t ia = 0; ia < NA; ia++)
      { if (seen[ia] == 0)
          { fprintf(wr, "    f%-7d =", kf);

            tosl_arc_id_t ka = ia;
            int32_t face_ok = 1;
            do 
              { /* Get and print vertex id: */
                tosl_vert_id_t kv = mesh->Arc[ka].ivorg;
                fprintf(wr, " v%d", kv);
                if ((kv < 0) || (kv >= mesh->NV)) 
                  { fprintf(wr, " !! invalid vertex");
                    face_ok = 0;
                    break;
                  }
                char *xka = tosl_arc_id_to_string(ka);
                fprintf(wr, " %s", xka);
                free(xka);
                if (seen[ka] != 0)
                  { fprintf(wr, " !! faces merge at arc %d", ka);
                    face_ok = 0;
                    break;
                  }

                seen[ka] = 1; 
                
                /* Advance to next arc: */
                tosl_arc_id_t ja = mesh->Arc[ka].skip;
                if (ja == -1) 
                  { fprintf(wr, " !! arc %d has no {.skip} link", ka);
                    face_ok = 0;
                    break;
                  }
                ka = ja;
                
              } while (ka != ia);
              
            if (face_ok != 0)
              { double area, nrm[3], ctr[3];
                tosl_mesh_face_normal_area_center(ia, mesh, &area, nrm, ctr);
                fprintf(wr, " area = %.2f normal = (", area);
                for (int32_t j = 0; j < 3; j++) { fprintf(wr, " %+6.4f", nrm[j]); }
                fprintf(wr, " ) center = (");
                for (int32_t j = 0; j < 3; j++) { fprintf(wr, " %+11.3f", ctr[j]); }
                fprintf(wr, " )");
              }
            fprintf(wr, "\n");
            kf++;
          }
      }
    fprintf(wr, "\n");
    fflush(wr);
    free(seen);
    return;
  }

void tosl_mesh_face_normal_area_center(tosl_arc_id_t ka, tosl_mesh_t *mesh, double *area_P, double nrm[], double ctr[])
  {
    int32_t debug = 0;
    
    double sum_nrm[3] = { 0.0, 0.0, 0.0 }; /* Sum of face normal times area of face. */
    double sum_ctr[3] = { 0.0, 0.0, 0.0 }; /* Sum of triangle barycenter times area. */

    auto void acc_normal_barycenter(tosl_vert_id_t iva, tosl_vert_id_t ivb, tosl_vert_id_t ivc);
      /* Gets the mesh vertices {va,vb,vc} with ids {iva,ivb,ivc}. Adds
        to {sum_nrm[0..2]} the area of triangle {va,vb,vc} times its
        normal direction. Namely, half the cross product of {vb-va} and
        {vc-va}. Also adds to {sum_ctr[0..2]} the coordinates of the
        triangle's barycenter times its area. */

    tosl_vert_id_t iv0 = -1; /* First vertex of face. */
    tosl_vert_id_t ivp = -1; /* Previous vertex of face other than {v0}. */
 
    tosl_arc_id_t ia = ka;
    tosl_arc_id_t ja = ka; /* Walks twice as fast around face, to catch bad loops. */
    int32_t face_ok = 1;
    int32_t deg = 0;
    while (1) 
      { /* Get the origin of {ia}: */
        tosl_vert_id_t iv = mesh->Arc[ia].ivorg;
        if ((iv < 0) || (iv >= mesh->NV)) 
          { if (debug) { fprintf(stderr, " !! invalid vertex index %d\n", iv); }  face_ok = 0; break; }
    
        /* Save the first and previous vertex: */
        if (iv0 == -1) 
          { iv0 = iv; } 
        else 
          { if (ivp != -1) 
              { /* Add area of triangle {v0,vp,vi}: */
                
                if (debug) { fprintf(stderr, "  accumulating v%d v%d v%d\n", iv0, ivp, iv); }
                acc_normal_barycenter(iv0, ivp, iv);
              }
            ivp = iv;
          }
        deg++;

        /* Advance {ia}: */
        ia = mesh->Arc[ia].skip;
        if (ia == -1) { if (debug) { fprintf(stderr, " !! missing {.skip} link\n"); } face_ok = 0; break; }

        if (ia == ka) { /* Completed face: */ break; }

        /* Advance {ja} by two steps: */
        ja = mesh->Arc[ja].skip;
        if (ja == -1) { if (debug) { fprintf(stderr, " !! missing {.skip} link\n"); } face_ok = 0; break; }
        ja = mesh->Arc[ja].skip;
        if (ja == -1) { if (debug) { fprintf(stderr, " !! missing {.skip} link\n"); } face_ok = 0; break; }
        if (ja == ia) { if (debug) { fprintf(stderr, " !! bizarre loop\n"); } face_ok = 0; break; }
      }
    if (debug) {  fprintf(stderr, "  degree = %d face_ok = %c\n", deg, "FT"[face_ok]); }

    if (face_ok != 0)
      { assert(ia == ka); /* Must have ended after a full round. */
        double sum_area2 = 0;
        for (int32_t j = 0; j < 3; j++) { sum_area2 += sum_nrm[j]*sum_nrm[j]; }
        double area = sqrt(sum_area2);
        for (int32_t j = 0; j < 3; j++) { nrm[j] = sum_nrm[j]/area; ctr[j] = sum_ctr[j]/area; }
        (*area_P) = area;
      }
    else
      { for (int32_t j = 0; j < 3; j++) { nrm[j] = NAN; ctr[j] = NAN; }
        (*area_P) = NAN;
      }
    return;

    void acc_normal_barycenter(tosl_vert_id_t iva, tosl_vert_id_t ivb, tosl_vert_id_t ivc)
      { tosl_point_t *va = &(mesh->Vpos[iva]);
        tosl_point_t *vb = &(mesh->Vpos[ivb]);
        tosl_point_t *vc = &(mesh->Vpos[ivc]);
        double u[3], v[3];
        for (int32_t j = 0; j < 3; j++)
          { u[j] = ((double)vb->c[j] - va->c[j]);
            v[j] = ((double)vc->c[j] - va->c[j]);
          }
        double cr[3]; /* Half of the cross product of {vb-va} and {vc-va}. */
        cr[0] = 0.5*(u[1]*v[2] - u[2]*v[1]);
        cr[1] = 0.5*(u[2]*v[0] - u[0]*v[2]);
        cr[2] = 0.5*(u[0]*v[1] - u[1]*v[0]);
        double area = sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]); /* Area of triangle. */
        for (int32_t j = 0; j < 3; j++)
          { sum_nrm[j] += cr[j];
            sum_ctr[j] += area*(va->c[j] + vb->c[j] + vc->c[j])/3;
          }
      }
  }

void tosl_mesh_coord_range_get(tosl_mesh_t *mesh, tosl_point_t *vmin_P, tosl_point_t *vmax_P)
  {
    tosl_point_t vmin = (tosl_point_t){{ INT32_MAX, INT32_MAX, INT32_MAX }};
    tosl_point_t vmax = (tosl_point_t){{ INT32_MIN, INT32_MIN, INT32_MIN }};
    for (tosl_vert_id_t iv = 0; iv < mesh->NV; iv++)
      { for (int32_t j = 0; j < 3; j++)
          { tosl_coord_t Cv = mesh->Vpos[iv].c[j];
            if (Cv > vmax.c[j]) { vmax.c[j] = Cv; }
            if (Cv < vmin.c[j]) { vmin.c[j] = Cv; }
          }
      }
    (*vmin_P) = vmin;
    (*vmax_P) = vmax;
  }
 
void tosl_mesh_coord_range_print(FILE *wr, char *pref, tosl_point_t *vmin, tosl_point_t *vmax, char *suff)
  {
    if (pref != NULL) { fputs(pref, wr); }
    for (int32_t j = 0; j < 3; j++)
      { fprintf(wr, "%s{%+d ..%+d}", (j == 0 ? "" : "×"), vmin->c[j], vmax->c[j]); }
    if (suff != NULL) { fputs(suff, wr); }
  }

tosl_mesh_t *tosl_mesh_new(int32_t NE_max, int32_t NV_max) 
  {
    int32_t NA_max = 2*NE_max;
    
    tosl_arc_t *Arc = malloc(NA_max*sizeof(tosl_arc_t));
    tosl_point_t *Vpos = malloc(NV_max*sizeof(tosl_point_t));
    char **Vlab = malloc(NV_max*sizeof(char*));

    tosl_mesh_t *mesh = malloc(sizeof(tosl_mesh_t));
    (*mesh) = (tosl_mesh_t) { .NE = 0, .Arc = Arc, .NV = 0, .Vpos = Vpos, .Vlab = Vlab, .NE_max = NE_max, .NV_max = NV_max };
    return mesh;
  }

void tosl_mesh_free(tosl_mesh_t *mesh)
  { if (mesh->Vlab != NULL)
      { for (tosl_vert_id_t iv = 0; iv < mesh->NV; iv++) { free(mesh->Vlab[iv]); }
        free(mesh->Vlab);
      }
    free(mesh->Vpos);
    free(mesh->Arc);
    free(mesh);
  }

tosl_vert_id_t tosl_mesh_add_vert
  ( tosl_point_t *v,
    char *lab,
    tosl_mesh_t *mesh
  )
  {
    tosl_vert_id_t kv = mesh->NV;
    assert((kv >= 0) && (kv < mesh->NV_max));
    mesh->Vpos[kv] = *v;
    if (mesh->Vlab != NULL) { mesh->Vlab[kv] = lab; } else { assert(lab == NULL); }
    mesh->NV = kv + 1;
    return kv;
  }

tosl_arc_id_t tosl_mesh_add_edge
  ( tosl_vert_id_t kv0,
    tosl_vert_id_t kv1,
    int32_t set_skip, 
    tosl_mesh_t *mesh
  )
  {
    uint8_t debug = 0;
    
    if (debug) { fprintf(stderr, "      adding edge from v%d to v%d\n", kv0, kv1); }
    assert((kv0 >= 0) && (kv0 < mesh->NV));
    assert((kv1 >= 0) && (kv1 < mesh->NV));
    assert(kv0 != kv1);
    
    int32_t ke = mesh->NE;
    assert((ke >= 0) && (ke < mesh->NE_max));
    tosl_arc_id_t ia = 2*ke;
    tosl_arc_id_t ja = ia+1;
    mesh->Arc[ia] = (tosl_arc_t){ .skip = -1, .ivorg = kv0, .pred = ia, .succ = ia };
    mesh->Arc[ja] = (tosl_arc_t){ .skip = -1, .ivorg = kv1, .pred = ia, .succ = ia };
    if (set_skip != 0) { mesh->Arc[ia].skip = ja; mesh->Arc[ja].skip = ia; }
    mesh->NE = ke + 1;
    if (debug) 
      { char *xia = tosl_arc_id_to_string(ia);
        char *xja = tosl_arc_id_to_string(ja);
        fprintf(stderr, "      added arcs %s and %s\n", xia, xja);
        free(xia); free(xja);
      }
    return ia;
  }
  
void tosl_mesh_splice(tosl_arc_id_t ia, tosl_arc_id_t ja, tosl_mesh_t *mesh)
  {
    int32_t NA = 2*mesh->NE;
    assert((ia >= 0) && (ia < NA));
    assert((ja >= 0) && (ja < NA));
    tosl_arc_id_t ia1 = mesh->Arc[ia].skip;
    tosl_arc_id_t ja1 = mesh->Arc[ja].skip;
    mesh->Arc[ia].skip = ja1;
    mesh->Arc[ja].skip = ia1;
  }

tosl_arc_id_t tosl_mesh_add_ring(int32_t n, char *pref, tosl_mesh_t *mesh)
  {
    assert(n >= 2); /* No loops please. */
    
    /* Create the vertices: */
    tosl_vert_id_t kv[n]; /* Indices of vertices. */
    tosl_point_t v = (tosl_point_t){{ 0, 0, 0 }}; /* Undefined coordinates. */
    for (int32_t i = 0; i < n; i++)
      { 
        char *lab = jsprintf("%s%d", pref, i);
        kv[i] = tosl_mesh_add_vert(&v, lab, mesh);
      }
    
    /* Create the edges: */
    tosl_arc_id_t ka[n]; /* Indices of arcs in one sense. */
    for (int32_t i = 0; i < n; i++)
      { tosl_vert_id_t kv0 = kv[i];
        tosl_vert_id_t kv1 = kv[(i+1) % n];
        ka[i] = tosl_mesh_add_edge(kv0, kv1, 0, mesh);
      }
      
    /* Connect the edges through the {.skip} links: */
    for (int32_t i = 0; i < n; i++)
      { tosl_arc_id_t ka0 = ka[i];
        tosl_arc_id_t ka1 = ka[(i+1) % n];
        mesh->Arc[ka0].skip = ka1;
        mesh->Arc[tosl_sym(ka1)].skip = tosl_sym(ka0);
      }
    return ka[0];   
  }
  
void tosl_mesh_add_path(int32_t n, tosl_arc_id_t ia0, tosl_arc_id_t ia1, char *pref, tosl_mesh_t *mesh)
  {
    assert(n >= 0);
    
    /* Edges whose origins are the bridgeheads: */
    tosl_arc_id_t ja0 = mesh->Arc[ia0].skip;
    tosl_arc_id_t ja1 = mesh->Arc[ia1].skip;

    /* Create the intermediate vertices: */
    tosl_vert_id_t kv[n]; /* Indices of vertices. */
    tosl_point_t v = (tosl_point_t){{ 0, 0, 0 }}; /* Undefined coordinates. */
    for (int32_t b = 0; b < n; b++)
      { 
        char *lab = jsprintf("%s%d", pref, b);
        kv[b] = tosl_mesh_add_vert(&v, lab, mesh);
      }
      
    /* Get initial and final verts of path: */
    tosl_vert_id_t iv0 = mesh->Arc[ja0].ivorg;
    tosl_vert_id_t iv1 = mesh->Arc[ja1].ivorg;
    
    /* Create the edges: */
    tosl_arc_id_t ka[n+1]; /* Indices of arcs in one sense. */
    tosl_vert_id_t kvp = iv0; /* Previous (origin) vertex. */
    for (int32_t b = 0; b <= n; b++)
      { tosl_vert_id_t kvb = (b == n ? iv1 : kv[b]); /* Next (dest) vertex. */
        ka[b] = tosl_mesh_add_edge(kvp, kvb, 0, mesh);
        kvp = kvb;
      }
      
    /* Connect the edges through the {.skip} links: */
    tosl_arc_id_t kap = ia0;
    tosl_arc_id_t rap = ja0;
    for (int32_t b = 0; b <= n; b++)
      { tosl_arc_id_t iab = ka[b];
        mesh->Arc[kap].skip = iab;
        mesh->Arc[tosl_sym(iab)].skip = rap;
        kap = iab;
        rap = tosl_sym(iab);
      }
     mesh->Arc[kap].skip = ja1;
     mesh->Arc[ia1].skip = rap;
   }

void tosl_mesh_link_triang(int32_t parity, tosl_arc_id_t ka0, tosl_arc_id_t ka1, tosl_arc_id_t ka2, tosl_mesh_t *mesh)
  { 
    int32_t debug = 1;
  
    if (debug) { tosl_tri_arc_id_print(stderr, "  linking triangle ", ka0, ka1, ka2, NULL); }
    if (debug) { fprintf(stderr, " parity %d", parity); }
    if (parity == 0)
      { tosl_arc_id_t ja0 = tosl_sym(ka0);
        tosl_arc_id_t ja1 = tosl_sym(ka1);
        tosl_arc_id_t ja2 = tosl_sym(ka2);
        ka0 = ja2;
        ka1 = ja1;
        ka2 = ja0;
        if (debug) { tosl_tri_arc_id_print(stderr, " = ", ka0, ka1, ka2, NULL); }
      }
    if (debug) { fprintf(stderr, "\n"); }
    if (debug) 
      { tosl_mesh_arc_print(stderr, "    ka0 = ", ka0, "\n", mesh);
        tosl_mesh_arc_print(stderr, "    ka1 = ", ka1, "\n", mesh);
        tosl_mesh_arc_print(stderr, "    ka2 = ", ka2, "\n", mesh);
      }
    /* The arcs should be still unlinked: */
    assert( mesh->Arc[ka0].skip == -1); mesh->Arc[ka0].skip = ka1;
    assert( mesh->Arc[ka1].skip == -1); mesh->Arc[ka1].skip = ka2;
    assert( mesh->Arc[ka2].skip == -1); mesh->Arc[ka2].skip = ka0;
  }


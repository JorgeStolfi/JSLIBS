/* Last edited on 2024-10-07 20:21:41 by stolfi */
/* Test of {tosl_build_lists_bin_sec.h} and {tosl_build_lists_hash.h} */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_pick_planes.h>
#include <tosl_build_lists_direct.h>
#include <tosl_build_lists_bin_sec.h>
#include <tosl_build_lists_hash.h>

#define N_METHODS 5
  /* Num of list building methods. The {method} parameter below ranges in  {0..N_METHODS-1}. */

#define N_TYPES 4
  /* Num of ypes of slicing plane sets. The {type} parameter below ranges in {0..N_TYPES-1}. */

int32_t main(int32_t argc, char **argv);
   
tosl_arc_id_t *do_build_lists(int32_t NP, tosl_coord_t Zplane[], tosl_mesh_t *mesh, int32_t method, int32_t debug);
  /* Builds the arc lists for the given {mesh} and {Zplane[0..NP-1]} using 
    the specified {method}, namely:
    
      Method 0: {tosl_build_lists_direct}.

      Method 1: {tosl_build_lists_bin_sec} with {use_bin=1}, {use_sec=0} (binary search only).
      
      Method 2: {tosl_build_lists_bin_sec} with {use_bin=0}, {use_sec=1} (secant search only).
      
      Method 3: {tosl_build_lists_bin_sec} with {use_bin=use_sec=1} (mixed binary+secant search).
      
      Method 4: {tosl_build_lists_hash}.
  */

void do_test_consistency(int32_t type, int32_t method, int32_t debug);
  /* Generates random meshes and sets of planes of the given {type}, builds lists 
    using the given {method}, and checks onsistency of the result. 
    If {debug} is 1 (true) uses few edges and planes and prints diagnostics. 
    If {debug} is 0 (false) uses lots of edges and planes and prints only a summary. */

void do_test_times(int32_t type);
  /* Measures the build times per edge with various methods on 
    plane sets of the given {type}.. */

tosl_mesh_t *make_pseudo_mesh(int32_t NE, int32_t NV, tosl_coord_t Zmin, tosl_coord_t Zmax);
  /* Creates a mesh structure just good enough for testing the list
    building functions. It has {NV} vertices with random even vertex
    coordinates between {Zmin} and {Zmax}, and {NE} edges connecting
    those vertices in pairs at random. The {.skip} links are all left
    undefined, but the {.pred} and {.succ} make each arc into a
    singleton circular list of arcs. */
  
void pick_vertices(tosl_mesh_t *mesh, tosl_coord_t Zmin, tosl_coord_t Zmax);
  /* Assumes that {mesh.NV} is zero.  Adds to it {NV = mesh.NV_max}
    vertices with random even coordinates. The {Z}-coordinates are
    all even, in random order, scattered between {Zmin} and {Zmax}. */

void pick_arcs(tosl_mesh_t *mesh);
  /* Assumes that {mesh.NE} is zero but {mesh.NV} is positive. Adds {NE
    = mesh.NE_max} edges that connect random pairs of distinct vertices.
    
    The {.skip} links are set to {-1}. The links {Arc[ia].pred} and
    {Arc[ia].succ} are set to {ia}, thus turning each arc into a
    singleton circular list. */
 
void check_list_correctness(tosl_arc_id_t L[], int32_t NP, tosl_coord_t Zplane[], tosl_mesh_t *mesh);
  /* Checks that the result {L[0..NP-1]} of {tosl_build_lists_bin_sec}
    or {tosl_build_lists_hash} with the arcs of the {mesh}, and slicing 
    planes {Zplane[0..NP-1]}. */

void get_build_time_per_edge
  ( tosl_mesh_t *mesh,
    int32_t NP,
    tosl_coord_t Zplane[],
    int32_t method,
    double *avg_P,
    double *dev_P
  );
  /* Returns the time of for building the edge lists with a specified {mehod},
    as per {do_build_lists}. */

void pick_plane_sets
  ( FILE *wr,
    int32_t NNP,
    int32_t NP_full,
    tosl_coord_t ZVmin,
    tosl_coord_t ZVmax,
    int32_t type,
    int32_t NP_set[],
    tosl_coord_t *Zplane_set[]
  );
  /* Creates {NNP} sets of slicing planes for timing tests. The {Z}
    coordinates of each set are in {Zmin .. Zmax} with a suitable gap
    between {Zmin} and the first plane, and between the last plane and
    {Zmax}.
    
    The {type} parameter determines the distribution of planes in 
    the interval. See {tosl_pick_planes} for the possible values and meanings.
    
    Stores the sets in {Zplane_set[0..NNP-1]} and the respective
    sizes in {NP_set[0..NNP-1]}.  Writes the statistics of each 
    plane set to {wr}. */
  
int32_t main(int32_t argc, char **argv)
  {
    srandom(4615);

    for (int32_t debug = 1; debug >= 0; debug--)
      { for (int32_t method = N_METHODS-1; method >= 0; method--)
          { for (int32_t type = 0; type < N_TYPES; type++)
             { fprintf(stderr, "\n");
               do_test_consistency(type, method, debug);
               fprintf(stderr, "\n");
             }
          }
      }
    
    for (int32_t type = 0; type < N_TYPES; type++) 
      { fprintf(stderr, "\n");
        do_test_times(type);
        fprintf(stderr, "\n");
      }

    return 0;
  }
     
void do_test_consistency(int32_t type, int32_t method, int32_t debug)
  {
    fprintf(stderr, "> --- %s --------------------\n", __FUNCTION__);
    fprintf(stderr, "  checking correctness with plane set type = %d method = %d debug = %c\n", type, method, "FT"[debug]);

    int32_t NE, NP; /* Number of edges, planes. */
    int32_t dZ; /* Approximate gap between planes. */
    if (debug)
      { /* Small problem, big gaps: */
        NE = 10; NP = 64; dZ = 1024; 
      }
    else
      { /* Big problem, small gaps: */
        NE = tosl_mesh_MAX_EDGES;
        NP = tosl_mesh_MAX_PLANES;
        dZ = 16;
      }
    int32_t NV = 2*NE; /* Max number of vertices. */
    fprintf(stderr, "  NE = %d NP = %d\n", NE, NP);

    /* Desired vertex {Z} range. Leave enough avg spacing for non-uniform {type}. */
    tosl_coord_t ZVmin = -dZ*NP;
    tosl_coord_t ZVmax = +dZ*NP;
    fprintf(stderr, "  nominal vertex {Z} range = {%+d .. %+d}\n", ZVmin, ZVmax);
    assert(ZVmax - ZVmin > 0); /* Check for overflow in {ZVmax-ZVmin}. */

    /* Create the vertices and arcs for the test: */
    tosl_mesh_t *mesh = make_pseudo_mesh(NE, NV, ZVmin, ZVmax);

    /* Slicing plane {Z} range. Leave vertices below and above all. */
    int32_t skosh = (ZVmax - ZVmin)/(NP+1); /* Gap to leave between {ZVmin,ZVmax} and plane {Z}s. */
    tosl_coord_t ZPmin = ZVmin + skosh; ZPmin = 2*((ZPmin - (ZPmin & 1))/2) + 1;
    tosl_coord_t ZPmax = ZVmax - skosh; ZPmax = 2*((ZPmax - (ZPmax & 1))/2) - 1;
    fprintf(stderr, "  skosh = %d nominal plane {Z} range = {%+d .. %+d}", skosh, ZPmin, ZPmax);
    fprintf(stderr, "  plane spacing type = %d\n", type);

    int32_t verbose = 1;
    tosl_coord_t *Zplane = tosl_pick_planes(NP, ZPmin, ZPmax, type, verbose);
    tosl_pick_planes_print_stats(stderr, NP, Zplane, ZVmin, ZVmax);

    if ((method == 0) && (type != 0) && (NP > 200))
      { fprintf(stderr, "    (skipped -- would take way too long)\n"); }
    else
      { tosl_arc_id_t *L = do_build_lists(NP, Zplane, mesh, method, debug); 
        check_list_correctness(L, NP, Zplane, mesh);
        free(L);
      }
    free(Zplane);
    tosl_mesh_free(mesh);
    fprintf(stderr, "< --- %s --------------------\n", __FUNCTION__);
  }

tosl_arc_id_t *do_build_lists(int32_t NP, tosl_coord_t Zplane[], tosl_mesh_t *mesh, int32_t method, int32_t debug)
  { tosl_arc_id_t *L = NULL;
    switch(method)
      { case 0: L = tosl_build_lists_direct (NP, Zplane, mesh, debug);       break; /* Direct indexing. */
        case 1: L = tosl_build_lists_bin_sec(NP, Zplane, mesh, 1, 0, debug); break; /* Binary only. */           
        case 2: L = tosl_build_lists_bin_sec(NP, Zplane, mesh, 0, 1, debug); break; /* Secant only. */           
        case 3: L = tosl_build_lists_bin_sec(NP, Zplane, mesh, 1, 1, debug); break; /* Mixed binary+secant. */   
        case 4: L = tosl_build_lists_hash   (NP, Zplane, mesh, debug);       break; /* Hash table. */            
        default: assert(0);
      }
    return L;
  }
   
void do_test_times(int32_t type)
  { 
    fprintf(stderr, "> --- %s --------------------\n", __FUNCTION__);

    char *fname = NULL;
    asprintf(&fname, "out/test_t%d.txt", type);
    FILE *wr = fopen(fname, "w");

    /* Generate a large mesh: */
    int32_t NE_full = tosl_mesh_MAX_EDGES/8;   /* Max number of edges. */
    int32_t NV_full = 2*NE_full;              /* Max number of vertices. */
    int32_t NP_full = tosl_mesh_MAX_PLANES/8;  /* Max number of planes. */
        
    /* Desired vertex {Z} range. Leave enough avg spacing for non-uniform {type}. */
    tosl_coord_t ZVmin = -128*NP_full;
    tosl_coord_t ZVmax = +128*NP_full;
    fprintf(wr, "  nominal vertex {Z} range = {%+d .. %+d}\n", ZVmin, ZVmax);
    assert(ZVmax - ZVmin > 0); /* Check for overflow in {ZVmax-ZVmin}. */

    /* Create the "mesh" (just a set of disconnected edges): */
    tosl_mesh_t *mesh = make_pseudo_mesh(NE_full, NV_full, ZVmin, ZVmax);
    
    /* Create test plane sets of various sizes: */
    int32_t NNP = 4; /* Number of test plane sets. */
    int32_t NP_set[NNP]; /* Size of each plane set. */
    tosl_coord_t *Zplane_set[NNP];
    pick_plane_sets(wr, NNP, NP_full, ZVmin, ZVmax, type, NP_set, Zplane_set);

    /* Timing tests for various subsets: */
    fprintf(wr, "Timing tests (μs/edge) plane set type = %d\n", type);
    fprintf(wr, "%10s %10s  %22s  %22s  %22s  %22s  %22s\n", "NE", "NP", "direct", "binary", "secant", "both", "hash");
    char *dash10 = "----------";
    char *dash22 = "----------------------";
    fprintf(wr, "%10s %10s  %22s  %22s  %22s  %22s  %22s\n", dash10, dash10, dash22, dash22, dash22, dash22, dash22);
    for (int32_t kp = 0; kp < NNP; kp++)
      { int32_t NP = NP_set[kp]; /* Num edges to use in this test: */
        tosl_coord_t *Zplane = Zplane_set[kp];
        int32_t NE = NE_full; /* Num edges to use in this test: */
        for (int32_t ke = 0; ke < 4; ke++)
          { mesh->NE = NE; /* Truncate the mesh to the first {NE} edges only. */
            fprintf(wr, "%10d %10d", NE, NP);
            for (int32_t method = 0; method < N_METHODS; method++)
              { if ((method == 0) && (type != 0) && (NP > 200))
                  { fprintf(wr, "  %22s", "(skipped)"); }
                else
                  { double avg, dev;
                    get_build_time_per_edge(mesh, NP, Zplane, method, &avg, &dev);
                    fprintf(wr, "  %12.3f ± %7.3f", avg, dev);
                  }
             }
            fprintf(wr, "\n");
            NE = NE / 4;
          }
        fprintf(wr, "\n");
      } 
    fclose(wr); free(fname);
    for (int32_t kp = 0; kp < NNP; kp++) { free(Zplane_set[kp]); }
    tosl_mesh_free(mesh);
    fprintf(stderr, "< --- %s --------------------\n", __FUNCTION__);
  }
        
void get_build_time_per_edge
  ( tosl_mesh_t *mesh,
    int32_t NP,
    tosl_coord_t Zplane[],
    int32_t method,
    double *avg_P,
    double *dev_P
  )
  {
    /* Truncate the edge list: */
    int32_t NT = 10;
    double time[NT];
    for (int32_t kt = 0; kt < NT; kt++)
      { double tstart = tosl_user_cpu_time_usec();
        tosl_arc_id_t *L = do_build_lists(NP, Zplane, mesh, method, 0);
        double tstop = tosl_user_cpu_time_usec();
        time[kt] = (tstop - tstart)/mesh->NE;
        free(L);
      }
    /* Compute average and deviation of times: */
    double avg, dev;
    tosl_compute_avg_dev(NT, time, &avg, &dev);
    (*avg_P) = avg;
    (*dev_P) = dev;
  }
    
void pick_plane_sets
  ( FILE *wr,
    int32_t NNP,
    int32_t NP_full,
    tosl_coord_t ZVmin,
    tosl_coord_t ZVmax,
    int32_t type,
    int32_t NP_set[],
    tosl_coord_t *Zplane_set[]
  )
  {
    int32_t NP = NP_full; /* Num edges to use in this test: */
    for (int32_t kp = 0; kp < NNP; kp++)
      { fprintf(wr, "  ----------------------------------------------------------------------\n");
        fprintf(wr, "  plane set %d\n", kp);
        /* Choose the {Z} range {ZPmin .. ZPmax} of the slicing planes: */
        /* Leave a fair fraction of the vertices below and above all planes. */
        int32_t skosh = (ZVmax - ZVmin)/(NP+1); /* Gap to leave between {ZVmin.ZVmax} and plane {Z}s. */
        tosl_coord_t ZPmin = ZVmin + skosh; ZPmin = 2*((ZPmin - (ZPmin & 1))/2) + 1;
        tosl_coord_t ZPmax = ZVmax - skosh; ZPmax = 2*((ZPmax - (ZPmax & 1))/2) - 1;

        NP_set[kp] = NP;
        int32_t verbose = 1;
        Zplane_set[kp] = tosl_pick_planes(NP, ZPmin, ZPmax, type, verbose);
        tosl_pick_planes_print_stats(wr, NP, Zplane_set[kp], ZVmin, ZVmax);
        NP = NP/4;
      }
    fprintf(wr, "  ----------------------------------------------------------------------\n");
  }

tosl_mesh_t *make_pseudo_mesh(int32_t NE, int32_t NV, tosl_coord_t ZVmin, tosl_coord_t ZVmax)
  { 
    tosl_mesh_t *mesh = tosl_mesh_new(NE, NV);
    pick_vertices(mesh, ZVmin, ZVmax);
    pick_arcs(mesh);
    return mesh;
  }

void pick_vertices(tosl_mesh_t *mesh, tosl_coord_t ZVmin, tosl_coord_t ZVmax)
  {
    int32_t NV = mesh->NV_max;
    for (tosl_vert_id_t iv = 0; iv < NV; iv++)
      { tosl_point_t p;
        for (int32_t j = 0; j < 3; j++)
          { double r = ((double)random())/((double)INT32_MAX);
            tosl_coord_t cj = ZVmin + (int32_t)(r*(ZVmax - ZVmin));
            /* Ensure coord is even: */
            cj = 2*(cj/2);
            p.c[j] = cj;
          }
        mesh->Vpos[iv] = p;
        mesh->Vlab[iv] = NULL;
      }
    mesh->NV = NV;
  }
 
void pick_arcs(tosl_mesh_t *mesh)
  { int32_t NE = mesh->NE_max;
    int32_t NV = mesh->NV;
    for (int32_t ke = 0; ke < NE; ke++)
      { tosl_arc_id_t ia = 2*ke;
        tosl_arc_id_t ja = ia + 1;
        /* Choose a pair of distinct random vertices as origins of {ia} and {ja}: */
        double ri = ((double)random())/((double)INT32_MAX);
        tosl_vert_id_t iv = (tosl_vert_id_t)floor(ri*NV);
        assert((iv >= 0) && (iv < NV));
        
        double rj = ((double)random())/((double)INT32_MAX);
        tosl_vert_id_t jv = (tosl_vert_id_t)floor(rj*(NV-1));
        if (jv >= iv) { jv++; }
        assert((jv >= 0) && (jv < NV));
        assert(iv != jv);
        
        mesh->Arc[ia] = (tosl_arc_t){ .skip = -1, .ivorg = iv, .pred = ia, .succ = ia};
        mesh->Arc[ja] = (tosl_arc_t){ .skip = -1, .ivorg = jv, .pred = ja, .succ = ja};
      }
    mesh->NE = NE;
  }
    
void check_list_correctness(tosl_arc_id_t L[], int32_t NP, tosl_coord_t Zplane[], tosl_mesh_t *mesh)
  {
    int32_t NA = 2*mesh->NE;
    
    /* Set all {.skip} linsk to {-1}: */
    for (tosl_arc_id_t ia = 0; ia < NA; ia++)
      { mesh->Arc[ia].skip = -1; }
      
    /* Check the lists {L[0..NP-1]} and set the {.skip} links to -2: */
    for (tosl_plane_id_t ip = 0; ip < NP; ip++)
      { tosl_coord_t Zp = Zplane[ip];
        /* Scan and check the list {L[ip]}: */
        tosl_arc_id_t ia = L[ip];
        if (ia != -1)
          { tosl_arc_id_t ka = ia;
            do 
              { tosl_coord_t Zorg = mesh->Vpos[mesh->Arc[ka].ivorg].c[2]; /* {Z} of origin of {ka}.  */
                tosl_arc_id_t ja = tosl_sym(ka);
                tosl_coord_t Zdst = mesh->Vpos[mesh->Arc[ja].ivorg].c[2]; /* {Z} of destination of {ka}.  */
                assert(Zorg < Zdst); /* Arc is upward. */
                assert((ip == 0) || (Zorg > Zplane[ip-1]));  /* Arc starts above {Zplane[ip-1]}. */
                assert(Zorg < Zp); /* Arc starts below {Zp}. */
                assert(Zdst > Zp); /* Arc ends above {Zp}. */
                assert(Zorg != Zp); /* No vertex-plane clashes. */
                assert(Zdst != Zp); /* No vertex-plane clashes. */
                mesh->Arc[ka].skip = -2;
                ka = mesh->Arc[ka].succ;
              }
            while (ka != ia);
          }
      }
      
    /* Check that all upgoing arcs that cross a plane are in some list: */
    for (tosl_arc_id_t ia = 0; ia < NA; ia++)
      { if (mesh->Arc[ia].skip == -1)
          { /* Arc {ia} is not in any {L[]} list. */
            tosl_coord_t Zorg = mesh->Vpos[mesh->Arc[ia].ivorg].c[2]; /* {Z} of origin of {ia}.  */
            tosl_arc_id_t ja = tosl_sym(ia);
            tosl_coord_t Zdst = mesh->Vpos[mesh->Arc[ja].ivorg].c[2]; /* {Z} of destination of {ia}.  */
            if ((Zorg < Zdst) && (Zorg < Zplane[NP-1]) && (Zdst > Zplane[0]))
              { /* Arc is upwards and its {Z}-range intersects the range {[Zplane[0] _ Zplane[NP-1]]}: */
                for (tosl_plane_id_t ip = 0; ip < NP; ip++)
                  { tosl_coord_t Zp = Zplane[ip];
                    assert((Zorg > Zp) || (Zdst < Zp)); /* Arc does not cross plane. */
                    assert(Zorg != Zp); /* No vertex-plane clashes. */
                    assert(Zdst != Zp); /* No vertex-plane clashes. */
                  }
              }
            else
              {  /* No vertex-plane clashes: */
                assert((Zorg >= Zdst) || (Zorg > Zplane[NP-1]) || (Zdst < Zplane[0])); 
              }
          }
      }
  }

/*  Last edited on 2024-10-09 13:59:45 by stolfi */
/* Test of {tosl_mesh_slice} in {totposlic.h} */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <tosl.h>
#include <tosl_slice.h>
#include <tosl_mesh.h>
#include <tosl_pick_planes.h>
#include <tosl_mesh_slice.h>
#include <tosl_mesh_make_keg.h>
#include <tosl_build_lists_hash.h>
 
#define M_kegs_sma 2
#define M_psets_sma 1
  /* Number of kegs and plane sets for small tests. */
  
#define M_kegs_big 12
#define M_psets_big 4
  /* Number of kegs and plane sets for big tests. */
  
#define NP_MAX 10000
  /* Max number of planes in any test. */
   
typedef struct keg_spec_t { int32_t NS; int32_t NR; int32_t NB; } keg_spec_t;
  /* Parameters that define a keg mesh. */

int32_t main(int32_t argc, char **argv);

void do_test_batch(char *tag, int32_t M_kegs, keg_spec_t keg_specs[], int32_t M_psets, int32_t NPs[], int32_t NT, int8_t debug);
  /* Calls {do_test} with the keg parameters {keg_specs[0..M_kegs-1]}, plane sets with sizes
    {NPs[0..M_psets-1]}, and {NT} repetitions of topological slice. */


void do_test
  ( FILE *wr, 
    int32_t NS,
    int32_t NR,
    int32_t NB,
    tosl_mesh_t *mesh,
    int32_t NP,
    tosl_coord_t ZH,
    int64_t NM_exp,
    int32_t NT,
    int8_t debug
  );
  /* Tests the procedure {tosl_mesh_slice} {NT} times, all with the 
    mesh created by {tosl_mesh_make_keg(NS,NR,NB,...)} but different sets
    sets of {NP} slicing planes strictly inside {-ZH .. +ZH}.  
    
    Prints to {wr} one line of the timing table with the test
    parameters and the average and deviation of the user CPU time spent
    slicing, excluding the time to build the initia arc lists and the
    slice post-processing time.
    
    The post-processing of each slice simply does {tosl_slice_free} on it. */

tosl_coord_t *pick_planes(int32_t NP, tosl_coord_t Zmin, tosl_coord_t Zmax, int32_t type, int8_t debug);
  /* Assumes that the vertex {Z}-range is {Zmin .. Zmaz}. Returns a
    vector of {NP} {Z}-coordinates of planes, all odd, strictly
    increasing, with spacing of type {type}, strictly inside the range
    {Zmin .. Zmax}. */

double check_time
  ( tosl_mesh_t *mesh,
    int32_t NP,
    tosl_coord_t Zplane[],
    tosl_arc_id_t L[] 
  );
  /* Returns the time of {tosl_mesh_slice}, in seconds, for the given
    {mesh} and the slicing planes with Z-coordinates {Zplane[0..NP-1]}.
    
    Assumes that {L[ip]}, for {ip} in {0..NP-1}, is the head arc id of
    the list of upward arcs tha start between {Zplane[ip-1]} and
    {Zplane[ip]} and end somewhere above the latter. */

void write_column_headers(FILE *wr);
  /* Writes to {wr} the headers of the timing table columns. */
  
int32_t main(int32_t argc, char **argv)
  {
    srandom(4615);
    
    keg_spec_t keg_specs_sma[M_kegs_sma] = {
        (keg_spec_t){ .NS =   4, .NR =  1, .NB =   0 },
        (keg_spec_t){ .NS =   4, .NR =  3, .NB =   2 },
      };
    int32_t NPs_sma[M_psets_sma] = { 7 };  /* Values of {NP} to use in the small tests. */

    keg_spec_t keg_specs_big[M_kegs_big] = {
        (keg_spec_t){ .NS =   10, .NR =  6, .NB =   2 },
        (keg_spec_t){ .NS =  250, .NR = 10, .NB =  50 },
        (keg_spec_t){ .NS =  250, .NR = 10, .NB = 100 },
        (keg_spec_t){ .NS =  250, .NR = 10, .NB = 200 },
        (keg_spec_t){ .NS =  250, .NR = 20, .NB =  50 },
        (keg_spec_t){ .NS =  250, .NR = 20, .NB = 100 },
        (keg_spec_t){ .NS =  250, .NR = 40, .NB =  50 },
        (keg_spec_t){ .NS =  500, .NR = 10, .NB =  50 },
        (keg_spec_t){ .NS =  500, .NR = 10, .NB = 100 },
        (keg_spec_t){ .NS =  500, .NR = 20, .NB =  50 },
        (keg_spec_t){ .NS = 1000, .NR = 10, .NB =  50 },
        (keg_spec_t){ .NS = 1500, .NR = 15, .NB =  75 }
      };
    int32_t NPs_big[M_psets_big] = { 1000, 2000, 4000, 8000 };  /* Values of {NP} to use in the big tests. */

    do_test_batch("sma", M_kegs_sma, keg_specs_sma, M_psets_sma, NPs_sma, 1, 1);
    do_test_batch("big", M_kegs_big, keg_specs_big, M_psets_big, NPs_big, 10, 0);
    
    return 0;
  }
  
void do_test_batch(char *tag, int32_t M_kegs, keg_spec_t keg_specs[], int32_t M_psets, int32_t NPs[], int32_t NT, int8_t debug)
  {
    char *fname = NULL; char *fname = jsprintf("out/times-%s.txt", tag);
    FILE *wr = fopen(fname, "w");
    free(fname); 

    write_column_headers(wr);
    for (uint32_t ikeg = 0;  ikeg < M_kegs; ikeg++)
      { /* Get the keg parameters: */
        keg_spec_t *spkeg = &(keg_specs[ikeg]);
        int32_t NS = spkeg->NS;
        int32_t NR = spkeg->NR;
        int32_t NB = spkeg->NB;
        
        /* Decide the total {Z} span of the keg: */
        int32_t dZ_min = 10; /* Min plane spacing and edge height. */
        tosl_coord_t Zspan_plane = dZ_min*(NP_MAX+1); /* So that the plane spacing is about {dZ_min}. */
        tosl_coord_t Zspan_edge = dZ_min*NR*(NB+1); /* So that the vertex spacing is about {dZ_min}. */
        tosl_coord_t Zspan = (Zspan_plane > Zspan_edge ? Zspan_plane : Zspan_edge);
        
        /* Create the test mesh: */
        tosl_coord_t ZH = Zspan/2; /* The keg will span from {-ZH} to {+ZH}. */
        fprintf(stderr, "creating the test mesh NS = %d  NR = %d  NB = %d  ZH = %d...\n", NS, NR, NB, ZH);
        tosl_mesh_t *mesh = tosl_mesh_make_keg(NS, NR, NB, ZH, debug);
        if (debug) { fprintf(stderr, "\ntest mesh:\n"); tosl_mesh_print(stderr, mesh);  fprintf(stderr, "\n\n"); }
        fprintf(stderr, "checking the mesh...\n");
        tosl_mesh_check(mesh); 

        /* Slice it with various plane set sizes: */
        for (uint32_t ipset = 0;  ipset < M_psets; ipset++)
          { int32_t NP = NPs[ipset];
            assert(NP < NP_MAX);

            /* Estimate number of plane-edge intersections: */
            int64_t NM_exp = ((int64_t)NP)*NS;
            fprintf(stderr, "  testing with NP = %d,  estimated NM = %ld\n", NP, NM_exp);

            do_test(wr, NS, NR, NB, mesh, NP, ZH, NM_exp, NT, debug);
          }
        fprintf(wr, "\n");
        fprintf(stderr, "\n\n======================================================================\n\n");
        tosl_mesh_free(mesh);
      }
    fclose(wr);
  }

void do_test
  ( FILE *wr,
    int32_t NS,
    int32_t NR,
    int32_t NB,
    tosl_mesh_t *mesh,
    int32_t NP,
    tosl_coord_t ZH,
    int64_t NM_exp,
    int32_t NT,
    int8_t debug
  )           
  {
    /* Timing tests: */
    double time[NT]; /* CPU times of runs. */
    for (uint32_t it = 0;  it < NT; it++)
      { fprintf(stderr, "choosing slicing planes...\n");
        int32_t type = 3; /* Inter-plane spacing distribution. */
        tosl_coord_t *Zplane = pick_planes(NP, -ZH, +ZH, type, debug);
        
        fprintf(stderr, "creating bucket lists...\n");
        tosl_arc_id_t *L = tosl_build_lists_hash(NP, Zplane, mesh, 0);
    
        /* Run the slicer: */
        int32_t NM_cmp = 0;  /* Total number of intersections actualy found. */
        
        auto void post(tosl_plane_id_t ip, tosl_slice_t *S);
          /* Consumes a slice.  Mainly, runs {free} on it. */
        
        fprintf(stderr, "slicing mesh...\n");
        tosl_mesh_slice(mesh, NP, Zplane, L, &post, &(time[it]), debug);
        
        assert(NM_cmp == NM_exp);
        free(Zplane);
        free(L);
        
        void post(tosl_plane_id_t ip, tosl_slice_t *S)
          { assert(S->Z == Zplane[ip]);
            NM_cmp += S->NV;
            tosl_slice_free(S);
          }
      }
    /* Compute average and deviation of times: */
    double avg_t, dev_t;
    tosl_compute_avg_dev(NT, time, &avg_t, &dev_t);

    char *fmt_data = "  %6d %6d %6d  %10d %10d  %10d  %12d  %8.3f %8.3f\n";
    
    fprintf(wr, fmt_data, NS, NR, NB, mesh->NE, mesh->NV, NP, NM_exp, avg_t/1.0e6, dev_t/1.0e6);

    write_column_headers(stderr);
    fprintf(stderr, fmt_data, NS, NR, NB, mesh->NE, mesh->NV, NP, NM_exp, avg_t/1.0e6, dev_t/1.0e6);
  }
            
void write_column_headers(FILE *wr)
  {
    fprintf(wr, "# Times for {tosl_slice} exclusive of list building (seconds)\n");
    char *dash6 = "------";
    char *dash8 = "--------";
    char *dash10 = "----------";
    char *dash12 = "------------";
    char *fmt_header = "# %6s %6s %6s  %10s %10s  %10s  %12s  %8s %8s\n";
    fprintf(wr, fmt_header, "NS", "NR", "NB", "NE", "NV", "NP", "NM", "avg", "dev");
    fprintf(wr, fmt_header, dash6, dash6, dash6, dash10, dash10, dash10, dash12, dash8, dash8);
  }    

tosl_coord_t *pick_planes(int32_t NP, tosl_coord_t Zmin, tosl_coord_t Zmax, int32_t type, int8_t debug)
  {
    int32_t skosh = (Zmax - Zmin)/(NP+1); /* Spacing between vertex {Z} extremes and plane {Z} extremes. */
    tosl_coord_t Z0 = Zmin + skosh;  Z0 = 2*((Z0 - (Z0&1))/2) + 1;
    tosl_coord_t Z1 = Zmax - skosh;  Z1 = 2*((Z1 - (Z1&1))/2) + 1;
    if (debug) { fprintf(stderr, "    Z0 = %+d  Z1 = %+d\n", Z0, Z1); }
    int32_t verbose = debug;
    tosl_coord_t *Zplane = tosl_pick_planes(NP, Z0, Z1, type, verbose);
    
    return Zplane;
  }

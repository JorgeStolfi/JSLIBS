#define PROG_NAME "test_maps_expand_reduce"
#define PROG_DESC "checks the {pst_{height,slope}_map_{expand,reduce}.h} routines"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-17 04:26:12 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_maps_expand_reduce_C_COPYRIGHT \
  "Copyright © 2010  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -function {ZFUNC} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  test_encode_gamma(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2010-08-14 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_maps_expand_reduce_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program checks the functions {pst_height_map_expand}, {pst_height_map_reduce}," \
  " {pst_slope_map_expand}, {pst_slope_map_reduce}."

#define PROG_INFO_OPTS \
  "  -function {ZFUNC}\n" \
  "    Specifies the number of the function to use. See {pst_proc_map.h}"

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <float_image.h>
#include <pst_slope_map.h>
#include <pst_height_map.h>
#include <pst_proc_map.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <vec.h>
#include <argparser.h>

typedef struct options_t
  { int32_t function;  /* Number of function to use. */
  } options_t;
  
typedef pst_proc_map_zfunc_t zfunc_t;

options_t *tb_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int32_t main(int32_t argc, char **argv);

void test_r_e_slope_map(int32_t nfunc);
  /* Tests {pst_slope_map_expand}, {pst_slope_map_reduce}.  The integer {nfunc} is 
    the iteration number. Output files will be called "out/{nfunc}-{NX}-{NY}-{tag}.fni"
    where {nfunc} is a procedural image ID number (see {pst_proc_map.h}),
    {NX,NY} are the dimensions of the test map, and {tag} is "IG", "JG", etc. */

void make_maps
  ( int32_t NX,
    int32_t NY,
    zfunc_t func,
    double maxGDiff,
    float_image_t **ZP,
    float_image_t **GP
  );
  /* Creates a height map {*ZP} of size {(NX+1)×(NY+1)} from a
    procedural image {func}. Also creates its gradient map {*GP}, both
    with size {NX×NY}. The {maxGDiff} is the threshold that defines the
    weight map (see {pst_proc_map_make_images}). Allocates the images
    internally. If any of {ZP,GP} is null, omits the corresponding
    map. */

void compare_and_write_slope_maps
  ( int32_t nfunc,
    float_image_t *IZ,
    float_image_t *IG, 
    float_image_t *JG, 
    float_image_t *KZ, 
    float_image_t *KG
  );
  /* Takes a gradient map {IG}, a reduced version {JG}
     of the same, and a reference image {KG} for {JG}. Writes all
     images and also the difference maps {DG=JG-KG}.  The parameter {nfunc} is the 
     function ID used for file names. All images are written to files 
     "out/{nfunc}-{NX}-{NY}-{tag}.fni" where {NX,NY} are the dimensions of {IG},
     and {tag} identifies the map. */

void write_map(float_image_t *M, int32_t nfunc, int32_t NX, int32_t NY, char *tag);
  /* Writes the image {M} as a FNI file called "out/{nfunc}-{NX}-{NY}-{tag}.fni", where
    {nfunc} is formatted as 4 digits zero-filled.  The parameters {NX} and {NY} 
    need not be the dimensions of {M}. */

void test_r_e_height_map(int32_t nfunc);
  /* Tests {pst_height_map_expand}, {pst_height_map_reduce}.  The integer {nfunc} is 
    the iteration number. Output files will be called "out/{nfunc}-{tag}.fni"
    where {tag} is "IZ", "IG", etc. */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tb_parse_options(argc, argv);
    
    /* Create a test map: */
    float_image_t *IZ = NULL;
    float_image_t *IG = NULL;
    test_r_e_slope_map(o->function);
    test_r_e_height_map(o->function);
    fprintf(stderr, "\n");
    float_image_free(IZ);
    float_image_free(IG);

    /* Cleanup: */
    free(o); o = NULL;
    return 0;
  }

void test_r_e_slope_map(int32_t nfunc)
  {
    fprintf(stderr, "--- testing slopes ---------------------------------------------------\n");
    /* Choose size: */
    int32_t NXI = int32_abrandom(10, 20);
    int32_t NYI = int32_abrandom(10, 20);
    /* Choose the function: */
    zfunc_t *func = pst_proc_map_function_generic(nfunc);  /* The function that defines the map. */
    double maxGDiff = NAN; /* Threshold for weight map definition. */
    switch(nfunc)
      { 
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
          maxGDiff = +INF;
          break;
        case 8:
          maxGDiff = 0.01;
          break;
        case 5:
        case 6:
        case 7:
        case 9:
        case 10:
        case 17:
        case 18:
        case 19:
        case 20:
          maxGDiff = 0.001;
          break;
        default:
          maxGDiff = +INF;
          break;
      }
    /* Create full-size height, slope and weight maps: */
    float_image_t *IZ = NULL;
    float_image_t *IG = NULL;
    make_maps(NXI, NYI, func, maxGDiff, &IZ, &IG);
    /* Shrink it: */
    float_image_t *JG = pst_slope_map_shrink(IG);
    /* Make a reduced version of the height, slope, and weight maps: */
    int32_t NX_J = (NXI + 1)/2;
    int32_t NY_J = (NYI + 1)/2;
    float_image_t *KZ = NULL;
    float_image_t *KG = NULL;
    make_maps(NX_J, NY_J, func, maxGDiff, &KZ, &KG);
    /* Compare and write: */
    compare_and_write_slope_maps(nfunc, IZ, IG, JG, KZ, KG);
    fprintf(stderr, "----------------------------------------------------------------------\n");
    float_image_free(JG);
    float_image_free(KG);
  }
  
void test_r_e_height_map(int32_t nfunc) { }

void compare_and_write_slope_maps
  ( int32_t nfunc,
    float_image_t *IZ, 
    float_image_t *IG, 
    float_image_t *JG, 
    float_image_t *KZ,
    float_image_t *KG
  )  
  {
    int32_t NCI, NXI, NYI;
    float_image_get_size(IG, &NCI, &NXI, &NYI);
    demand(NCI == 3, "gradient map has wrong depth");
    
    /* Dimensions of reduced reference map: */
    int32_t NCK, NXK, NYK;
    float_image_get_size(IG, &NCK, &NXK, &NYK);
    demand(NCK == 3, "reduced gradient map has wrong depth");
    
    /* Use the full dimensions to compose the filename: */
    write_map(IZ, nfunc, NXI, NYI, "IZ");
    write_map(IG, nfunc, NXI, NYI, "IG");
    write_map(JG, nfunc, NXI, NYI, "JG");
    write_map(KZ, nfunc, NXI, NYI, "KZ");
    write_map(KG, nfunc, NXI, NYI, "KG");
    
    /* Compare the reduced gradient maps (shrunk and pristine): */
    float_image_t *DG = float_image_new(3, NXK, NYK);
    float_image_mix_channels(1.0,JG,0, -1.0,KG,0, DG,0);
    float_image_mix_channels(1.0,JG,1, -1.0,KG,1, DG,1);
    float_image_mix_channels(1.0,JG,2, -1.0,KG,2, DG,2);
    write_map(DG, nfunc, NXI, NYI, "DG");
    
    float_image_free(DG);
  }
  
void write_map(float_image_t *M, int32_t nfunc, int32_t NX,int32_t NY, char *tag)
  { 
    char *fname = jsprintf("out/%02d-%04d-%04d-%s.fni", nfunc, NX, NY, tag);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, M);
    fclose(wr);
    free(fname);
  }

void make_maps
  ( int32_t NX,
    int32_t NY,
    zfunc_t func,
    double maxGDiff,
    float_image_t **ZP,
    float_image_t **GP
  )
  {
    int32_t smoothGL = 2;
    int32_t smoothGN = 11;
    pst_proc_map_sampling_t smpG = pst_proc_map_make_sampling_tables(smoothGL, smoothGN);
    int32_t smoothZL = 1;
    int32_t smoothZN = 11;
    pst_proc_map_sampling_t smpZ = pst_proc_map_make_sampling_tables(smoothZL, smoothZN);
    bool_t numGrad = FALSE;
    double sigmaG = 0.0;
    double sigmaW = 0.0;
    float_image_t *Z = float_image_new(2, NX+1, NY+1); /* Must be non-null. */
    float_image_t *G = (GP == NULL ? NULL : float_image_new(3, NX, NY));
    pst_proc_map_make_images(func, NX, NY, smpZ, smpG, numGrad, maxGDiff, sigmaG, sigmaW, Z, G, NULL);
    if (ZP != NULL) { (*ZP) = Z; }
    if (GP != NULL) { (*GP) = G; }
    pst_proc_map_free_sampling_tables(&smpZ);
    pst_proc_map_free_sampling_tables(&smpG);
  }

options_t *tb_parse_options(int32_t argc, char **argv)
  {
    /* INITIALIZATION: */

    /* Start the command line analyzer {pp}: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);

    /* Process "-help" and "-info" options: */
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    options_t *o = (options_t *)malloc(sizeof(options_t));

    /* PARSE KEYWORD ARGUMENTS: */
    argparser_get_keyword(pp, "-function");
    o->function = (int32_t)argparser_get_next_int(pp, pst_proc_map_MIN_ZFUNC, pst_proc_map_MAX_ZFUNC);
    
    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }

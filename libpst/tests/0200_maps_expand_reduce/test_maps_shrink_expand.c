#define PROG_NAME "test_maps_shrink_expand"
#define PROG_DESC "checks the {pst_{vertex,cell}_map_{expand,shrink}.h} routines"
#define PROG_VERS "1.0"

/* Last edited on 2025-03-01 21:02:48 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_maps_shrink_expand_C_COPYRIGHT \
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
  "  " test_maps_shrink_expand_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program checks the functions {pst_vertex_map_expand}, {pst_vertex_map_shrink}," \
  " {pst_cell_map_expand}, {pst_cell_map_shrink}.\n" \
  "\n" \
  "  Output files will be called \"{outPrefix}-{mapTag}{dirTag}{varTag}.fni\" where\n" \
  "\n" \
  "    {outPrefix} is \"out/{func_num:02d}-{func_name}-{NXI:04d}-{NYI:04d}\"\n" \
  "\n" \
  "    {func_num} and {func_name} are the height function's number and name (see {pst_proc_map.h})\n" \
  "\n" \
  "    {NXI,NYI} are the dimensions of the test map before shrinking or expanding\n" \
  "\n" \
  "    {mapTag} is \"GI\" for slope maps and \"ZI\" for height maps\n" \
  "\n" \
  "    {dirTag} is \"\" for the original map, \"S\" for the shrunk version, \"E\" for" \
  " the expanded version\n" \
  "\n" \
  "    {varTag} is "" for the original map, \"-cmp\" for the result of shrink or" \
  " expand, \"-ref\" for the ideal version of the same, and \"-dif\" for the" \
  " difference between \"*-cmp\" and \"*-ref\""

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
#include <pst_height_map.h>
#include <pst_slope_map.h>
#include <pst_cell_map_shrink.h>
#include <pst_vertex_map_shrink.h>
#include <pst_vertex_map_expand.h>
#include <pst_proc_map.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <vec.h>
#include <wt_table_hann.h>
#include <argparser.h>

typedef struct options_t
  { int32_t function;  /* Number of function to use. */
  } options_t;
  
typedef pst_proc_map_zfunc_t zfunc_t;

options_t *tb_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int32_t main(int32_t argc, char **argv);

void test_pst_cell_map_shrink(int32_t nfunc);
  /* Tests {pst_cell_map_shrink}, {pst_cell_map_expand}.  The integer {nfunc} is 
    the function number. */

void test_pst_vertex_map_shrink__pst_vertex_map_expand(int32_t nfunc);
  /* Tests {pst_vertex_map_shrink} and {pst_vertex_map_expand}. The integer
    {nfunc} is the function number. */

void tmse_do_test_vertex_map_shrink_expand
  ( float_image_t *ZI,
    pst_proc_map_zfunc_props_t *fp,
    r2_t *org,
    double pixSize,
    bool_t shrink
  );
  /* Tests {pst_vertex_map_shrink} or {pst_vertex_map_expand} on the
    image {ZI} depending on {shrink}. Compares the shrunk/expanded height
    map {ZR_cmp} with an height maps {ZR_exp} of the same size created
    from scratch with {fp}. The {org} and {pixSize}
    are the parameters used to build {ZI}. */
    
float_image_t* tmse_make_cell_map
  ( int32_t NXG,
    int32_t NYG,
    r2_t *org,
    double pixSize,
    pst_proc_map_zfunc_props_t *fp,
    int32_t xDebug,
    int32_t yDebug
  );
  /* Creates a gradient map {GI} with 3 channel and size {NXG×NYG}.
    The {maxGrad} and {maxGDiff} parameters are passed on to
    {pst_proc_map_make_slope_map}. */

float_image_t* tmse_make_vertex_map
  ( int32_t NXG,
    int32_t NYG,
    r2_t *org,
    double pixSize,
    pst_proc_map_zfunc_props_t *fp,
    int32_t xDebug,
    int32_t yDebug
  );
  /* Creates a height map {ZI} with 2 channels and size {(NXG+1)×(NYG+1)} from a
    procedural image {func}, using {pst_proc_map_make_height_map}. */

void tmse_compare_and_write_maps
  ( int32_t func_num,
    char *func_name,
    int32_t NCI,
    int32_t NXI,
    int32_t NYI,
    char *mapTag,
    float_image_t *M_cmp, 
    float_image_t *M_ref
  );
  /* Takes a map {M_cmp} that was obtained by reducing or
    expanding some original map {M} with {NCI} channels and size {NXI × NYI}, and a
    reference image {M_ref} for {M_cmp}. The {dirTag} should be "S" for
    shrink and "E" for expand.  Writes these two images images and also
    the difference maps {M_dif = M_cmp - M_ref}. The parameter {nfunc} is
    the function ID used for file names. Mll images are written to files
    "out/{nfunc:02d}-{NXI:04d}-{NYI:04d}-{mapTag}-{verTag}.fni" and {verTag} is
    "cmp", "ref", or "dif". */

void tmse_write_map(float_image_t *M, int32_t func_num, char *func_name, int32_t NCI, int32_t NXI, int32_t NYI, char *tag, char *sub);
  /* Writes the image {M} as a FNI file called "out/{nfunc:02d}-{NXI:04d}-{NYI:04d}-{tag}-[sub}.fni", where
    {nfunc} is formatted as 4 digits zero-filled.  The image must have {NCI} channels,
    but its dimensions need not be {NXI} and {NYI}. The "-{sub}" part is omitted if {sub} is {NULL}. */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tb_parse_options(argc, argv);
    
    /* Create a test map: */
    test_pst_cell_map_shrink(o->function);
    test_pst_vertex_map_shrink__pst_vertex_map_expand(o->function);

    /* Cleanup: */
    free(o); o = NULL;
    return 0;
  }

void test_pst_cell_map_shrink(int32_t nfunc)
  {
    fprintf(stderr, "--- testing {pst_cell_map_shrink} ---------------------------------------------------\n");
    /* Choose slope map size: */
    int32_t NC = 3;
    int32_t NXGI = 2*int32_abrandom(5, 10);
    int32_t NYGI = 2*int32_abrandom(5, 10)+1;

    /* Choose the function: */
    pst_proc_map_zfunc_props_t fp = pst_proc_map_function_generic(nfunc);  /* The function that defines the map. */

    /* Create full-size slope map: */
    int32_t NGMIN = (NXGI < NYGI ? NXGI : NYGI);  /* Smallest dimension of slope map. */
    double pixSize = 2.0/NGMIN; /* Number of grid pixels per {func} slope map domain unit. */
    r2_t org = (r2_t){{ 0.5*NXGI, 0.5*NYGI }};
    float_image_t *GI = tmse_make_cell_map(NXGI, NYGI, &org, pixSize, &fp, -1, -1);
    tmse_write_map(GI, fp.num, fp.name, NC, NXGI, NYGI, "GI", NULL);

    /* Shrink it: */
    int32_t NXGS = (NXGI + 1)/2; /* Expected shrunk cols. */
    int32_t NYGS = (NYGI + 1)/2; /* Expected shrunk rows. */
    float_image_t *GS_cmp = pst_cell_map_shrink(GI, 2, NXGS, NYGS, 1.0);

    /* Make a reduced version of the slope map from scratch: */
    double scaleXY = 0.5;
    r2_t org_ref; r2_scale(scaleXY, &org, &org_ref);
    double pixSize_ref = pixSize*scaleXY;
    float_image_t *GS_ref = tmse_make_cell_map(NXGS, NYGS, &org_ref, pixSize_ref, &fp, -1, -1);

    /* Compare and write: */
    tmse_compare_and_write_maps(fp.num, fp.name, NC, NXGI, NYGI, "GS", GS_cmp, GS_ref);

    float_image_free(GS_cmp);
    float_image_free(GS_ref);
    fprintf(stderr, "----------------------------------------------------------------------\n");
  }
   
void test_pst_vertex_map_shrink__pst_vertex_map_expand(int32_t nfunc) 
  { 
    fprintf(stderr, "--- testing {pst_vertex_map_shrink} ---------------------------------------------------\n");

    /* Choose slope map size: */
    int32_t NC = 2;
    int32_t NXGI = 2*int32_abrandom(5, 10);
    int32_t NYGI = 2*int32_abrandom(5, 10)+1;

    /* Choose the function: */
    pst_proc_map_zfunc_props_t fp = pst_proc_map_function_generic(nfunc);  /* The function that defines the map. */

    /* Create full-size height map: */
    int32_t NGMIN = (NXGI < NYGI ? NXGI : NYGI);  /* Smallest dimension of slope map. */
    double pixSize = 2.0/NGMIN; /* Number of grid pixels per {func} slope map domain unit. */
    r2_t org = (r2_t){{ 0.5*NXGI, 0.5*NYGI }};
    float_image_t *ZI = tmse_make_vertex_map(NXGI, NYGI, &org, pixSize, &fp, 4, 4);
    tmse_write_map(ZI, fp.num, fp.name, NC, NXGI, NYGI, "ZI", NULL);

    tmse_do_test_vertex_map_shrink_expand(ZI, &fp, &org, pixSize, TRUE);
    tmse_do_test_vertex_map_shrink_expand(ZI, &fp, &org, pixSize, FALSE);
    fprintf(stderr, "----------------------------------------------------------------------\n");
  }
    
void tmse_do_test_vertex_map_shrink_expand
  ( float_image_t *ZI,
    pst_proc_map_zfunc_props_t *fp,
    r2_t *org,
    double pixSize,
    bool_t shrink
  ) 
  { 
    /* Get size of initial height map {ZI}: */
    int32_t NC, NXZI, NYZI;
    float_image_get_size(ZI, &NC, &NXZI, &NYZI);
    demand((NC == 1) || (NC == 2), "bad depth of height map {ZI}");
    
    /* Dimensions of the corresponding slope map: */
    int32_t NXGI = NXZI-1;
    int32_t NYGI = NYZI-1;
    
    /* Shrink or expand {ZI}: */
    int32_t NXZR, NYZR; /* Expected cols and rows of shrunk/expanded height map {ZR}. */
    float_image_t *ZR_cmp = NULL; /* Shrunk/expanded height map. */
    char *tag = NULL; /* "ZS" or "ZE". */
    double scale = NAN; /* Scaling factor for {Z} and {XY}. */
    int32_t xDebug, yDebug;
    if (shrink)
      { NXZR = NXZI/2 + 1;
        NYZR = NYZI/2 + 1;
        fprintf(stderr, "shrinking from %d×%d to %d×%d ...\n", NXZI, NYZI, NXZR, NYZR);
        scale = 0.5;
        ZR_cmp = pst_vertex_map_shrink(ZI, 1, NXZR, NYZR, scale);
        tag = "ZS";
        xDebug = 2; yDebug = 2;
      }
    else
      { NXZR = 2*NXZI-1 + 0; 
        NYZR = 2*NYZI-1 + 1; 
        fprintf(stderr, "expanding from %d×%d to %d×%d ...\n", NXZI, NYZI, NXZR, NYZR);
        scale = 2.0;
        ZR_cmp = pst_vertex_map_expand(ZI, 1, NXZR, NYZR, scale);
        tag = "ZE";
        xDebug = 8; yDebug = 8;
      }
    float_image_check_size(ZR_cmp, NC, NXZR, NYZR, "shrunk/expanded heghth map has wrong dept or size");

    /* Make a srunk/expanded version of the height map from scratch: */
    r2_t org_ref; r2_scale(scale, org, &org_ref);
    double pixSize_ref = pixSize*scale;
    float_image_t *ZR_ref = tmse_make_vertex_map(NXZR-1, NYZR-1, &org_ref, pixSize_ref, fp, xDebug, yDebug);

    /* Compare and write: */
    tmse_compare_and_write_maps(fp->num, fp->name, NC, NXGI, NYGI, tag, ZR_cmp, ZR_ref);

    float_image_free(ZR_cmp);
    float_image_free(ZR_ref);
  }

void tmse_compare_and_write_maps
  ( int32_t func_num,
    char *func_name,
    int32_t NC, 
    int32_t NXGI, 
    int32_t NYGI,
    char *tag,
    float_image_t *M_cmp, 
    float_image_t *M_ref
  )  
  {
    /* Dimensions of computed map: */
    int32_t NCM, NXM, NYM;
    float_image_get_size(M_cmp, &NCM, &NXM, &NYM);
    demand(NCM == NC, "computed map {M_cmp} map has wrong depth");
    
    /* Dimensions of reference map: */
    float_image_check_size(M_ref, NCM, NXM, NYM, "reference map {M_ref} has wrong depth or size");
    
    /* Use the dimensions of original cell map to compose the filename: */
    tmse_write_map(M_cmp, func_num, func_name, NC, NXGI, NYGI, tag, "cmp");
    tmse_write_map(M_ref, func_num, func_name, NC, NXGI, NYGI, tag, "ref");
    
    /* Compare the reduced cell maps (computed and expected): */
    float_image_t *M_dif = float_image_new(NC, NXM, NYM);
    for (int32_t c = 0; c < NCM; c++)
      { float_image_mix_channels(1.0,M_cmp,c, -1.0,M_ref,c, M_dif,c); }
    tmse_write_map(M_dif, func_num, func_name, NC, NXGI, NYGI, tag, "dif");
    
    float_image_free(M_dif);
  }

void tmse_write_map(float_image_t *M, int32_t func_num, char *func_name, int32_t NC, int32_t NX, int32_t NY, char *tag, char *sub)
  { int32_t NCM, NXM, NYM;
    float_image_get_size(M, &NCM, &NXM, &NYM);
    demand(NCM == NC, "map {M} has wrong depth");
    char *fullTag = (sub == NULL ? tag : jsprintf("%s-%s", tag, sub));
    char *fname = jsprintf("out/%02d-%s-%04d-%04d-%s.fni", func_num, func_name, NX, NY, fullTag);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, M);
    fclose(wr);
    free(fname);
    if (fullTag != tag) { free(fullTag); }
  }

float_image_t* tmse_make_cell_map
  ( int32_t NXG,
    int32_t NYG,
    r2_t *org,
    double pixSize,
    pst_proc_map_zfunc_props_t *fp,
    int32_t xDebug,
    int32_t yDebug
  )
  {
    uint32_t NS = 5;
    double ws[NS];
    uint32_t stride;
    wt_table_hann_fill(NS, 0.0, ws, &stride);
    assert(stride == NS/2 + 1);

    bool_t numGrad = TRUE;
    double minWeight = 1.0e-4;
    float_image_t *G = pst_proc_map_make_slope_map
      ( fp->func, NXG, NYG, org, pixSize, NS, ws, numGrad, fp->maxGrad, fp->maxGDiff, minWeight, xDebug, yDebug );
    float_image_check_size(G, 3, NXG, NYG, "{pst_proc_map_make_slope_map} returns wrong depth/size");
    return G;
  }

float_image_t *tmse_make_vertex_map
  ( int32_t NXG,
    int32_t NYG,
    r2_t *org,
    double pixSize,
    pst_proc_map_zfunc_props_t *fp,
    int32_t xDebug,
    int32_t yDebug
  )
  {
    uint32_t NS = 5;
    double ws[NS];
    uint32_t stride;
    wt_table_hann_fill(NS, 0.0, ws, &stride);
    assert(stride == NS/2 + 1);

    float_image_t *Z = pst_proc_map_make_height_map
      ( fp->func, NXG, NYG, org, pixSize, NS, ws, xDebug, yDebug );
    float_image_check_size(Z, 2, NXG+1, NYG+1, "{pst_proc_map_make_height_map} returns wrong depth/size");
    return Z;
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

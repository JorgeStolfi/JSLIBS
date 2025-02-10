#define PROG_NAME "test_mfok_render"
#define PROG_DESC "test of {multifok_rendere.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-02-08 17:28:11 by stolfi */
/* Created on 2025-02-02 by J. Stolfi, UNICAMP */

#define test_mfok_render_COPYRIGHT \
  "Copyright © 2025  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>

#include <bool.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <r2.h>
#include <r3.h>
#include <i2.h>
#include <affirm.h>
#include <argparser.h>
#include <float_image.h>
#include <float_image_test.h>
#include <float_image_interpolate.h>
#include <float_image_read_gen.h>
#include <float_image_write_gen.h>

#include <multifok_sampling.h>
#include <multifok_raytrace.h>
#include <multifok_frame.h>
#include <multifok_render.h>
#include <multifok_scene_object_raytrace.h>

#define PROG_INFO_OFILES \
  "Each run of the program generates a bunch of test image files.  Each" \
  " test image shows a bunch of balls rendered with a specific amount" \
  " of isotropic (\"ambient\") lighting and specific surface gloss" \
  " and lambedo (Lambertian albedo) patterns.  Each" \
  " ball is viewed from a point at infinity in the {Z=(0,0,1)} direction" \
  " and is illuminated by a white infinitely distant point-like source plus some" \
  " isotropic light, with the same color and intensity at every surface point" \
  " of each ball.  The direction of the point source makes an angle" \
  " with the vertical that varies from ball to ball, from zero (straight" \
  " down from above) to 180 degrees (straight up from below).\n" \
  "\n" \
  "   Each test image is written to a separate" \
  " file \"{frameFolder}/image.txt\".  " mren_write_frame_folder_INFO "\n" \
  "\n" \
  "   " mren_make_and_write_frame_code_INFO ""

typedef struct mren_options_t
  { int32_t imageSize_X;    /* Image width. */
    int32_t imageSize_Y;    /* Image height. */
    uint32_t pixSampling;   /* Pixel subsampling order. */
    i2_t debugPixel;        /* Pixel to debug. */
    char *outFolder;        /* Directory where to put output frame folders. */
  } mren_options_t;
  /* Command line parameters. */

int32_t main(int32_t argn, char **argv);

mren_options_t *mren_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

void mren_make_and_write_frame
  ( int32_t NX,
    int32_t NY,
    int32_t M,
    char LIT,
    char ISO,
    char GLO,
    char LAM,
    float_image_t *sGlo_pat,
    float_image_t *sLam_pat,
    uint32_t HS,
    char *outFolder,
    i2_t *iPixDeb
  );
  /* Creates an RGB test image {oimg} with size {M*NX,M*NY}, normally showing an {M×M}
    array of balls rendered with various lighting and surface finish parameters
    determinde by {LIT}, {ISO}, {GLO}, and {LAM} as described in 
    {mren_make_and_write_frame_code_INFO}.  
    
    The viewing direction is {(0,0,-1)} for all the balls. The direction
    {uLit} of the light source lies in some fixed vertical plane and varies
    from ball to ball, from {(0,0,+1)} to {(0,0,-1)} in increments of
    {180/24} degrees.
    
    If {LIT} is '0', the direction of the light source is not relevant.
    Therefore, in that case the frame shows a single ball, instead of an
    {M×M} array of balls.
    
    Patterned images (code 'P') use the pattern image {pat}.

    The image is written to file "{frameFolder}/{image}.png". See
    {mren_write_frame_folder_INFO} for the  meaning of {frameFolder}. 
    
    Uses {(2*HS+1)^2} subsampling points at each
    pixel. If {HS} is zero, each pixel of {oimg} will have
    only one sampling point, at its center. Otherwise the sampling
    points will extend outside that pixel, with 2D Hann window
    weighting. */
    
#define mren_write_frame_folder_INFO \
  "The folder name {frameFolder} is \"{outFolder}/{sizeTag}/{renderTag}\"," \
  " where {sizeTag} is \"nb{MM}-{XXXX}x{YYYY}\", {renderTag}" \
  " is \"/lit{LIT}-iso{ISO}-glo{GLO}-lam{LAM}\", " \
  " {MM} is the value of {M} formatted as '%02d', {XXXX} and {YYYY} are the image dimensions" \
  " as '%04d', {LIT} and {ISO} indicate the amount and color of the" \
  " point source and isotropic (\"ambient\") light, respectively; {GLO} indicates the" \
  " amount and color of the surface gloss; and {LAM} indicates the" \
  " underlying lambedo (Lambertian albedo)."
    
#define mren_make_and_write_frame_code_INFO \
  "Each of {LIT}, {ISO}, {GLO}, and {LAM} may" \
  " be '0', '1', 'C', or 'P'.  The codes '0' and '1' mean that the" \
  " respective parameter (point source intensity, isotropic light" \
  " instensity, gloss, or lambedo) is" \
  " uniform {(0,0,0)} or {(1,1,1)} over the whole sphere.  The code 'C' means that the" \
  " parameter is a random triplet between those two extremes, the" \
  " same at every point of the surface of every sphere.  The code 'P' (not" \
  " allowed for the light parameters {LIT} and {ISO}) means" \
  " that the parameter varies between those two extremes that varies over" \
  " the surface, but is the same at the same point of each sphere."

void mren_paint_frame_sub_image
  ( multifok_frame_t *fr,
    int32_t xLo, 
    int32_t xHi, 
    int32_t yLo,
    int32_t yHi,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
    multifok_pattern_frgb_proc_t *get_sGlo,
    multifok_pattern_frgb_proc_t *get_sLam,
    multifok_sampling_t *samp,
    multifok_raytrace_debug_pred_t *debug_pixel
  );
  /* Paints a sub-image of frame {fr} spanning cols {xLo..xHi}
    and rows {yLo..yHi}.  Assumes that the directional and ambient
    light colors are {cLit} and {cIso}, and that the 
    glossy and lambertian BRDF coeffs at a generic scene point
    {p} are given by {get_sGlo(p)} and {get_sLam(p)}, respectively. */

multifok_frame_t*  mren_make_frame(int32_t M, int32_t NX, int32_t NY);
  /* Creates a frame {fr} where every image consists of
    {M×M} sub-images of {NX} cols and {NY} rows.
    Will have {fr.zFoc = 0} and {fr.zDep = +INF} */

void mren_write_frame(multifok_frame_t *fr, char *frameFolder);
  /* Writes selected images from frame {fr} to files
    "{frameFolder}/{imageName}.png" where {imageName} is "sVal", "sNrm", "hAvg", etc. */

void mren_write_image(float_image_t *oimg, float vMin, float vMax, char *frameFolder, char *imageName);
  /* Writes the image {oimg} to file  "{frameFolder}/{imageName}.png",
    mapping {[vMin _ vMax]} to {0..maxval}. . */
    
float_image_t *mren_read_pattern_image(char *patName);
  /* Reads an image from file "in/{patName}.png". */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    fprintf(stderr, "starting %s ...\n", PROG_NAME);
    
    mren_options_t *o = mren_parse_options(argc, argv);
    
    int32_t NX = o->imageSize_X;
    int32_t NY = o->imageSize_Y;
    
    int32_t M = 5;
    uint32_t HS = o->pixSampling;
    
    /* Read the pattern image: */
    float_image_t *sGlo_pat = mren_read_pattern_image("NOISE06");
    float_image_t *sLam_pat = mren_read_pattern_image("melon26");

    i2_t *iPixDeb = &(o->debugPixel);

    /* Colored ambient light plus colored unidir light, patterned Lambertian albedo and glossy component: */
    mren_make_and_write_frame(NX, NY, M, 'C', 'C', 'P', 'P', sGlo_pat, sLam_pat, HS, o->outFolder, iPixDeb); 

    /* White unidir light, no ambient light, no Lambertian albedo, white glossy component: */
    /* mren_make_and_write_frame(NX, NY, M, '1', '0', '1', '0', sGlo_pat, sLam_pat, HS, o->outFolder, iPixDeb);  */

    /* White ambient light, no unidir light, no Lambertian albedo, white glossy component: */
    /* mren_make_and_write_frame(NX, NY, M, '0', '1', '1', '0', sGlo_pat, sLam_pat, HS, o->outFolder, iPixDeb); */ 

    /* White unidir light, no ambient light, no glossy component, white Lambertian albedo: */
    /* mren_make_and_write_frame(NX, NY, M, '1', '0', '0', '1', sGlo_pat, sLam_pat, HS, o->outFolder, iPixDeb);  */

    /* White ambient light, no unidir light, no glossy component, white Lambertian albedo: */
    /* mren_make_and_write_frame(NX, NY, M, '0', '1', '0', '1', sGlo_pat, sLam_pat, HS, o->outFolder, iPixDeb);  */

    /* Colored ambient light plus colored unidir light, no glossy component, patterned Lambertian albedo: */
    /* mren_make_and_write_frame(NX, NY, M, 'C', 'C', '0', 'P', sGlo_pat, sLam_pat, HS, o->outFolder, iPixDeb);  */
    
    return 0;
  }
  
void mren_make_and_write_frame
  ( int32_t NX,
    int32_t NY,
    int32_t M,
    char LIT,
    char ISO,
    char GLO,
    char LAM,
    float_image_t *sGlo_pat,
    float_image_t *sLam_pat,
    uint32_t HS,
    char *outFolder,
    i2_t *iPixDeb
  )
  {
    fprintf(stderr, "generating image pixels\n");

    multifok_sampling_t *samp = multifok_sampling_choose(HS, 1, FALSE);
    /* Paranoia: */
    assert(samp->NS >= 1);
    assert(samp->KR == 1);
   
    double azimLit = M_PI/3; /* Azimuth of point light source. */
    
    auto double eval_spherical_pattern(r3_t *p, float_image_t *pat);
      /* A "procedural" image function that actually samples the monochrome image {pat} at 
        the scene point {p} with bilinear interpolation. */
 
    auto bool_t debug_pixel(i2_t *iPix);
      /* A procedure suitable for the {debug_pixel} argument of
        {multifok_raytrace_paint_frame_rectangle}. It returns {TRUE} iff
        the pixel {iPix} is congruent to {iPixDeb} modulo {NX,NY}. */
    
    auto frgb_t get_param(r3_t *q, char code, frgb_t *fg, frgb_t *bg, float_image_t *pat, double rot);
      /* Returns either {(0,0,0)}, {(1,1,1)}, {fg}, or a pattern-directed 
        mix of {bg} and {fg}, depending on {code} being '0', '1', 'C', or 'P',
        respectively.  If {code} is not 'P', the parameters {q,bg,fpat,rot} are
        ignored and may be {NULL}. */

    demand(LIT != 'P', "invalid {LIT} code");
    frgb_t cLit_fg = (frgb_t) {{ 0.900f, 0.600f, 0.100f }};
    frgb_t cLit = get_param(NULL, LIT, &cLit_fg, NULL, NULL, 0); /* Unidirectional light color. */
    
    demand(ISO != 'P', "invalid {ISO} code");
    frgb_t cIso_fg = (frgb_t) {{ 0.100f, 0.300f, 0.100f }};
    frgb_t cIso = get_param(NULL, ISO, &cIso_fg, NULL, NULL, 0); /* Isotropic light color. */

    frgb_t sGlo_fg = (frgb_t) {{ 0.300f, 0.700f, 0.900f }}; /* For when {GLO} is 'C' or 'P'. */
    frgb_t sGlo_bg = (frgb_t) {{ 0.050f, 0.030f, 0.000f }}; /* For when {GLO} is 'P'. */
    
    /* If {LAM} is 'C' or 'P': */
    frgb_t sLam_fg = (frgb_t) {{ 0.900f, 0.500f, 0.200f }}; /* For when {LAM} is 'C' or 'P'. */
    frgb_t sLam_bg = (frgb_t) {{ 0.300f, 0.250f, 0.050f }}; /* For when {LAM} is 'P'. */

    auto frgb_t get_sGlo(r3_t *q);
    auto frgb_t get_sLam(r3_t *q);
      /* Obtain the {sGlo} and {sLam} surface finish parameters given
        the respective codes {GLO,LAM} and the position of point {(x,y,z)}
        relative to the center, normalized to unit length. */
    
    multifok_frame_t *fr = mren_make_frame(M, NX, NY);
    
    if (frgb_is_all_zeros(&cLit))
      { /* No unidirectional light, so there is no point in varying {uLit}: */
        int32_t xLo = 0, xHi = M*NX - 1;
        int32_t yLo = 0, yHi = M*NY - 1;
        r3_t uLit = (r3_t){{ 0,0,0 }}; /* Shoul dnot matter. */
        mren_paint_frame_sub_image
          ( fr, xLo, xHi, yLo, yHi,
            &uLit, &cLit, &cIso, &get_sGlo, &get_sLam,
            samp, debug_pixel
          );
      }
    else
      { /* Vary {uLit} from dead front to dead back: */
        uint32_t kxy = 0; /* Sub-images generated so far. */
        r3_t uLit; /* Point light source direction for current sub-image. */
        for (int32_t ky = 0; ky < M; ky++)
          { for (int32_t kx = 0; kx < M; kx++)
              { /* Render sub-image in column {kx} and row {ky}: */
                int32_t xLo = kx*NX, xHi = xLo + NX - 1;
                int32_t yLo = ky*NY, yHi = yLo + NY - 1;

                double elevLit = M_PI/2*(1.0 - (2.0*kxy)/(M*M));

                uLit = (r3_t){{ cos(azimLit)*cos(elevLit), sin(azimLit)*cos(elevLit), sin(elevLit) }};

                mren_paint_frame_sub_image
                  ( fr, xLo, xHi, yLo, yHi,
                    &uLit, &cLit, &cIso, &get_sGlo, &get_sLam,
                    samp, debug_pixel
                  );

                kxy++;
              }
           }
       }

    char *sizeTag = jsprintf("nb%02d-%04dx%04d", M, NX, NY);
    char *sizeFolder = jsprintf("%s/%s", outFolder, sizeTag);

    char *renderTag = jsprintf("lit%c-iso%c-glo%c-lam%c", LIT, ISO, GLO, LAM);
    char *frameFolder = jsprintf("%s/%s", sizeFolder, renderTag);

    mkdir(sizeFolder, 0755);
    mkdir(frameFolder, 0755);

    mren_write_frame(fr, frameFolder);
    
    free(sizeTag); free(renderTag);
    free(sizeFolder); free(frameFolder);
    multifok_sampling_free(samp);
    multifok_frame_free(fr);
    return;
    
    double eval_spherical_pattern(r3_t *p, float_image_t *pat)
      { assert(pat != NULL);
        /* Get pattern image dimensions: */
        int32_t NCP, NXP, NYP;
        float_image_get_size(pat, &NCP, &NXP, &NYP);
        demand(NCP == 1, "patterns must be monochrome");
        demand(NXP == NYP, "patterns must be square");
        int32_t NXYP = NXP;
        /* Project {p} onto the unit sphere: */
        r3_t u3; (void)r3_dir(p, &u3);
        /* Mirror the bottom hemisphere to the upper one: */ 
        u3.c[2] = fabs(u3.c[2]);
        /* Stereographic projection of {u} to the {Z=0} plane from the {(0,0,-1)} pole: */
        double s = 1 + u3.c[2];
        r2_t u2 = (r2_t) {{ u3.c[0]/s, u3.c[1]/s }}; 
        /* Map {u2} from {[-1 _ +1]^2} to {[0 _ NXYP]}: */
        u2.c[0] = 0.5*(u2.c[0] + 1.0)*NXYP;
        u2.c[1] = 0.5*(u2.c[1] + 1.0)*NXYP;
        /* Fetch the pattern image, phew: */
        int32_t order = 1;
        double res = float_image_interpolate_sample(pat, 0, u2.c[0], u2.c[1], order, ix_reduce_mode_PXMIRR);
        return res;
      }

    frgb_t get_sGlo(r3_t *q)
      { return get_param(q, GLO, &sGlo_fg, &sGlo_bg, sGlo_pat, 0); }
      
    frgb_t get_sLam(r3_t *q)
      { return get_param(q, LAM, &sLam_fg, &sLam_bg, sLam_pat, M_PI/3); }
      
    frgb_t get_param(r3_t *q, char code, frgb_t *fg, frgb_t *bg, float_image_t *pat, double rot)
      { 
        if (code == '0')
          { return (frgb_t){{ 0,0,0 }}; }
        else if (code == '1')
          { return (frgb_t){{ 1,1,1 }}; }
        else if (code == 'C')
          { return *fg; }
        else if (code == 'P')
          { r3_t u = *q; 
            if (rot != 0)
              { r3_rot_axis(&u, 0, 1, rot, &u);
                r3_rot_axis(&u, 1, 2, rot, &u);
                r3_rot_axis(&u, 2, 0, rot, &u);
              }
            double r = eval_spherical_pattern(&u, pat);
            return frgb_mix(1-r, fg, r, bg);
          }
        else
          { demand(FALSE, "invalid render code"); }
      }

    bool_t debug_pixel(i2_t *iPix)
      { bool_t debug_col = (iPix->c[0] - iPixDeb->c[0]) % NX == 0;
        bool_t debug_row = (iPix->c[1] - iPixDeb->c[1]) % NY == 0;
        bool_t debug = (debug_col && debug_row);
        return debug;
      }
  }

void mren_paint_frame_sub_image
  ( multifok_frame_t *fr,
    int32_t xLo, 
    int32_t xHi, 
    int32_t yLo,
    int32_t yHi,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
    multifok_pattern_frgb_proc_t *get_sGlo,
    multifok_pattern_frgb_proc_t *get_sLam,
    multifok_sampling_t *samp,
    multifok_raytrace_debug_pred_t *debug_pixel
  )
  {
    fprintf(stderr, "generating sub-image [%d..%d, %d..%d]\n", xLo, xHi, yLo, yHi);
    
    /* Get frame dimensions: */
    int32_t NCF, NXF, NYF;
    float_image_get_size(fr->sVal, &NCF, &NXF, &NYF);
    demand((0 <= xLo) && (xLo < xHi) && (xHi < NXF), "invalid col range");
    demand((0 <= yLo) && (yLo < yHi) && (yHi < NYF), "invalid row range");
    double xCtr = 0.5*(xLo + xHi + 1);
    double yCtr = 0.5*(yLo + yHi + 1);

    r3_t dRef = (r3_t){{ 0,0,-1 }};
    
    int32_t NC = NCF; /* Channels in {sVal} image. */
    int32_t NX = xHi - xLo + 1;
    int32_t NY = yHi - yLo + 1;
    double NXY_min = (NX < NY ? NX : NY); 
    double xy_scale = 0.95*NXY_min; /* Image to scene scaling factor. */

    auto void trace_ray(r3_t *pR, r3_t *dR, bool_t debug_loc, r3_t *pHit_P, r3_t *sNrm_P, int32_t NC_loc, float sVal[]);
      /* A ray-tracing func suitable for {multifok_raytrace_paint_frame_rectangle}.
        The procedure assumes that the projection of the scene (the unit ball of {\RR^3})
        on the image plane is the 
        circle of radius {xy_scale} centered in the rectangle
        {[xLo _ xHi+1] × [yLo _ yHi+1]}. */
    
    auto void map_point(r2_t *p2_img, r3_t *p3_scene);
      /* An image-to-scene coordinate conversion function, suitable for
        the {map} argument of {multifok_raytrace_paint_frame_rectangle}. The image {X} and {Y}
        coordinates of {p2_img} are mapped to scene coordinates described above,
        and the {Z} corodinate is set to {zFoc}. */
    
    double zDep = fr->zDep;
    double zFoc = fr->zFoc;
    
    /* Bounding box of the ball: */
    interval_t bbox[3]; 
    for (int32_t j = 0; j < 3; j++) { bbox[j] = (interval_t){{ -1, +1 }}; }

    bool_t verbose_pixel = FALSE;
    multifok_raytrace_paint_frame_rectangle
      ( fr, xLo, xHi, yLo, yHi, trace_ray, map_point,
        &dRef, zFoc, zDep, samp,
        verbose_pixel, debug_pixel, NULL, NULL
      );

    return;

    void trace_ray(r3_t *pR, r3_t *dR, bool_t debug_loc, r3_t *pHit_P, r3_t *sNrm_P, int32_t NC_loc, float sVal[])
      { assert(NC_loc == NC);
        assert(dR->c[2] < 0);
        
        double debug = FALSE; /* !!! Fix this !!! */
        double tHit = multifok_scene_object_raytrace_BALL(bbox, pR, dR, debug); 
        r3_t pHit, sNrm;
        if (isfinite(tHit))
          { r3_mix(1.0, pR, tHit, dR, &pHit);
            /* The hit point must be on the upper half of the unit sphere: */
            demand(fabs(r3_norm(&pHit) - 1.0) < 1.0e-7, "ball tracing failed");
            demand(pHit.c[2] >= bbox[2].end[0]-1.0e-10, "ball tracing failed");
            pHit.c[2] = fmax(0, pHit.c[2]);
            /* The normal is just the hit point: */
            sNrm = pHit;
            /* Get the local gloss and lambedo: */
            frgb_t sGlo = get_sGlo(&pHit);
            frgb_t sLam = get_sLam(&pHit);
            /* Render the point: */
            frgb_t clr = multifok_render_compute_color
              ( dR, &sNrm, &sGlo, &sLam, uLit, cLit, cIso, debug ); 
            for (int32_t j = 0; j < NC; j++) { sVal[j] = (j < 3 ? clr.c[j] : clr.c[2]); }
          }
        else
          { r3_mix(1.0, pR, 3.0, dR, &pHit);
            sNrm = (r3_t){{ 0,0,0 }};
            for (int32_t j = 0; j < NC; j++) { sVal[j] = 0.500; }
          }

        (*pHit_P) = pHit;
        (*sNrm_P) = sNrm;
      }

    void map_point(r2_t *p2_img, r3_t *p3_scene)
      {
        double x_scene = 2*(p2_img->c[0] - xCtr)/xy_scale;
        double y_scene = 2*(p2_img->c[1] - yCtr)/xy_scale;
        (*p3_scene) = (r3_t){{ x_scene, y_scene, zFoc }};
      }
  }
    
multifok_frame_t*  mren_make_frame(int32_t M, int32_t NX, int32_t NY)
  { 
    int32_t NXF = M*NX;
    int32_t NYF = M*NY;
    int32_t NC = 3;

    double zDep = +INF;
    double zFoc = 0.0;
    
    float_image_t *sVal = float_image_new(NC, NXF, NYF);
    float_image_t *shrp = float_image_new(1, NXF, NYF);
    float_image_t *hAvg = float_image_new(1, NXF, NYF);
    float_image_t *hDev = float_image_new(1, NXF, NYF);
    float_image_t *sNrm = float_image_new(3, NXF, NYF);
  
    multifok_frame_t *fr = multifok_frame_from_images
      ( sVal,shrp,hAvg,hDev,sNrm, zFoc,zDep );
    return fr;
  }

void mren_write_frame(multifok_frame_t *fr, char *frameFolder)
  { mren_write_image(fr->sVal, 00.000f, +1.000f, frameFolder, "sVal");
    mren_write_image(fr->sNrm, -1.000f, +1.000f, frameFolder, "sNrm");
    mren_write_image(fr->hAvg, -1.001f, +1.001f, frameFolder, "hAvg");
  }

void mren_write_image(float_image_t *oimg, float vMin, float vMax, char *frameFolder, char *imgName)
  { char *fname = jsprintf("%s/%s.png", frameFolder, imgName);
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = TRUE;
    double gammaEnc = 1.0;
    double bias = 0.0;
    bool_t verbose = TRUE;
    float_image_write_gen_named(fname, oimg, ffmt, yUp, vMin, vMax, gammaEnc,bias, verbose);
    free(fname);
  }
  
float_image_t *mren_read_pattern_image(char *patName)
  { 
    char *fileName = jsprintf("in/texture/%s.png", patName);
    bool_t yUp = FALSE;
    double gammaDec, bias;
    bool_t verbose = TRUE;
    float_image_t *img = float_image_read_gen_named
      ( fileName, image_file_format_PNG, yUp, 
        0.0, 1.0, NULL, &gammaDec, &bias, 
        verbose
      );
    free(fileName);
    return img;
  }

mren_options_t *mren_parse_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);

    mren_options_t *o = notnull(malloc(sizeof(mren_options_t)), "no mem");

    argparser_get_keyword(pp, "-imageSize");
    o->imageSize_X = (int32_t)argparser_get_next_int(pp, 30, 4096);
    o->imageSize_Y = (int32_t)argparser_get_next_int(pp, 30, 4096);

    argparser_get_keyword(pp, "-pixSampling");
    o->pixSampling = (uint32_t)argparser_get_next_int(pp, 0, 9999);

    argparser_get_keyword(pp, "-debugPixel");
    o->debugPixel.c[0] = (int32_t)argparser_get_next_int(pp, 0, o->imageSize_X-1);
    o->debugPixel.c[1] = (int32_t)argparser_get_next_int(pp, 0, o->imageSize_Y-1);

    argparser_get_keyword(pp, "-outFolder");
    o->outFolder = argparser_get_next_non_keyword(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }

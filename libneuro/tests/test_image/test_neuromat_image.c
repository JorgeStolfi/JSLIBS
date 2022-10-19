#define PROG_NAME "test_neuromat_image"
#define PROG_DESC "Tests {neuromat_image.h}."
#define PROG_VERS "2021-08-24"

#define test_neuromat_image_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"
/* Last edited on 2021-08-29 01:07:31 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <float_image.h>
#include <frgb.h>
#include <affirm.h>
#include <jsfile.h>

#include <neuromat_image.h>
#include <neuromat_image_png.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_io.h>
 
float_image_t *tni_read_image(char *fname);
  /* Reads the FNI image file "fname". */
 
void tni_read_dataset(char *fname, double ***valP, int32_t *ntP, int32_t *ncP, int32_t *neP, int32_t *nmP);

int32_t main(int32_t argc, char **argv)
  {
    /* Read the scalp mask image {msk}, grab main image dimensions {NXI,NYI}: */
    float_image_t *msk = tni_read_image("in/280_320_e128_msk.fni");
    assert(msk->sz[0] == 1);
    int32_t NX = (int32_t)msk->sz[1];
    int32_t NY = (int32_t)msk->sz[2];
    
    /* Read a signed overlay image: */
    float_image_t *ovr = tni_read_image("in/280_320_e128_elp.fni");
    assert(ovr->sz[0] == 1);
    assert(ovr->sz[1] == NX);
    assert(ovr->sz[2] == NY);
 
    /* Read a potential field image */
    float_image_t *fld = tni_read_image("in/280_320_e128_fld.fni");
    assert(fld->sz[0] == 1);
    assert(fld->sz[1] == NX);
    assert(fld->sz[2] == NY);
    
    /* Create the ouptut image: */
    float_image_t *cim = float_image_new(4,NX,NY);
    double gamma = NAN; /* Use dafault gamma for output images. */
    
    fprintf(stderr, "--- testing {neuromat_image_colorize_field} ---\n");
    float_image_fill(cim, NAN);
    float vlo = +INF, vhi = -INF, vmax;
    float_image_update_sample_range(fld, 0, &vlo, &vhi);
    vmax = fmaxf(fabsf(vlo),fabsf(vhi));
    int32_t style = 1;
    neuromat_image_colorize_field(cim, fld, msk, (double)vmax, style);
    neuromat_image_png_write("out", "czf", cim, 0.0, 1.0, gamma);
    
    fprintf(stderr, "--- testing {neuromat_image_colorize_signed_overlay} ---\n");
    float_image_fill(cim, NAN);
    vlo = +INF; vhi = -INF;
    float_image_update_sample_range(ovr, 0, &vlo, &vhi);
    vmax = fmaxf(fabsf(vlo),fabsf(vhi));
    demand(vmax <= 1.0, "invalid signed overlay image");
    for (int32_t side = -1; side <= +1; side += 2)
      { frgb_t rgb = (side == 0 ? (frgb_t){{ 0.0, 0.5, 1.0 }} : (frgb_t){{ 1.0, 0.5, 0.0 }});
        neuromat_image_colorize_signed_overlay(cim, ovr, (sign_t)(2*side-1), rgb);
        char *tag = (side == 0 ? "cov_neg" : "cov_pos");
        neuromat_image_png_write("out", tag, cim, 0.0, 1.0, gamma);
      }
 
    fprintf(stderr, "--- {testing neuromat_image_paint_slider}  ---\n");
    float_image_fill(cim, 0.000f);
    vlo = +INF; vhi = -INF;
    int32_t xlo = 20;
    int32_t xsz = NX - 2*xlo;
    int32_t y0 = 20;
    int32_t y1 = y0 + 20;
    int32_t y2 = y1 + 20;
    frgb_t ftrack = (frgb_t){{ 0.500f, 0.600f, 0.700f }};
    frgb_t fslid = (frgb_t){{ 1.000f, 0.000f, 0.500f }};
    int32_t hw = 2;
    neuromat_image_paint_time_track(cim, hw, xlo, xsz, y0, ftrack);
    neuromat_image_paint_time_track(cim, hw, xlo, xsz, y1, ftrack);
    neuromat_image_paint_time_track(cim, hw, xlo, xsz, y2, ftrack);
    double tlo = 0.0;
    double thi = 500.0;
    double t = 0.67*tlo + 0.33*thi;
    int32_t ylo = y0;
    int32_t yhi = y2;
    int32_t shh = 4*hw;
    neuromat_image_paint_slider(cim, hw, shh, tlo, t, thi, xlo, xsz, ylo, yhi, fslid);
    neuromat_image_png_write("out", "sld", cim, 0.0, 1.0, gamma);

    float_image_free(msk);
    float_image_free(ovr);
    float_image_free(cim);
    
    return 0;
  }

float_image_t *tni_read_image(char *fname)
  { FILE *rd = open_read(fname, TRUE);
    float_image_t *msk = float_image_read(rd);
    fclose(rd);
    return msk;
  }
 
void tni_read_dataset(char *fname, double ***valP, int32_t *ntP, int32_t *ncP, int32_t *neP, int32_t *nmP)
  { 
    FILE *rd = open_read(fname, TRUE);
    
    /* Read the EEG header: */
    int32_t nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(rd, 20, 600, &nl);
    int32_t nc = h->nc;    /* Number of channels per data frame (including markers etc.). */
    int32_t ne = h->ne;    /* Number of electrodes. */
    int32_t nm = nc - ne;  /* Number of marker channel dots. */
    fprintf(stderr, "input file has %d channels\n", nc);
    fprintf(stderr, "channels comprise %d electrodes and %d marker channels\n", ne, nm);
    
    /* Read the EEG data frames: */
    /* We must read the frames with step 1, otherwise we may miss trigger pulses. */
    /* But we assume that the electrodes are filtered so that downsampling is aliasing-free.  */
    int32_t nt = 0;
    int32_t skip = 0;
    int32_t read = 0;
    double **val = neuromat_eeg_data_read(rd, skip, read, nc, &nl, &nt);
    fprintf(stderr, "read %d lines, got %d data frames\n", nl, nt);
    fclose(rd);
    
    (*ntP) = nt;
    (*ncP) = nc;
    (*neP) = ne;
    (*nmP) = nm;
    (*valP) = val;
  }


#define PROG_NAME "test_neuromat_image"
#define PROG_DESC "Tests {neuromat_image.h}."
#define PROG_VERS "2021-08-24"

#define test_neuromat_image_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"
/* Last edited on 2021-08-25 02:48:41 by stolfi */

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
    
    /* Read an EEG dataset: */
    double **val = NULL;
    int32_t nt, nc, ne, nm;
    tni_read_dataset("in/s019_r00418.txt", &val, &nt, &nc, &ne, &nm);
    assert(nm >= 2);
    
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
    neuromat_image_png_write("out/t", "czf", cim, 0.0, 1.0, gamma);
    
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
        neuromat_image_png_write("out/t", tag, cim, 0.0, 1.0, gamma);
      }
 
    fprintf(stderr, "--- testing painting time tracks and markers ---\n");
    float_image_fill(cim, 0.000f);
    vlo = +INF; vhi = -INF;
    int32_t xlo = 20;
    int32_t xsz = NX - 2*xlo;
    int32_t ylo = 20;
    int32_t y = ylo;
    frgb_t ftrack = (frgb_t){{ 0.500f, 0.600f, 0.700f }};
    frgb_t fmark = (frgb_t){{ 1.000f, 0.800f, 0.000f }};
    frgb_t fslid = (frgb_t){{ 1.000f, 0.000f, 0.500f }};
    int32_t hw = 2;
    int32_t ic_mark = nc - 1;
    neuromat_image_paint_time_track(cim, hw, xlo, xsz, y, ftrack);
    neuromat_image_paint_marker_ranges(cim, hw, nt, nc, val, ic_mark, xlo, xsz, y, fmark);
    double tlo = 0.0;
    double thi = (double)nt;
    double t = 0.67*tlo + 0.33*thi;
    int32_t yhi = ylo + 20;
    neuromat_image_paint_slider(cim, hw, tlo, t, thi, xlo, xsz, ylo, yhi, fslid);
    neuromat_image_png_write("out/t", "sld", cim, 0.0, 1.0, gamma);

    float_image_free(msk);
    float_image_free(ovr);
    for (int32_t it = 0; it < nt; it++) { free(val[it]; }
    free(val);
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


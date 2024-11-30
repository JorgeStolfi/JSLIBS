#define PROG_NAME "test_neuromat_eeg_image"
#define PROG_DESC "Tests {neuromat_eeg_image.h}."
#define PROG_VERS "2021-08-24"

#define test_neuromat_eeg_image_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-10-22 03:23:16 by stolfi */

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
#include <neuromat_eeg_image.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_io.h>
 
void tnei_read_dataset(char *fname, double ***valP, int32_t *ntP, int32_t *ncP, int32_t *neP, int32_t *nmP, char ***chnameP);
float_image_t **tnei_read_basis_images(char *prefix, int32_t ne, char *chname[], int32_t NX, int32_t NY);
void tnei_write_fni(char *tag, float_image_t *img);

int32_t main(int32_t argc, char **argv)
  {
    
    /* Read an EEG dataset: */
    double **val = NULL;
    int32_t nt, nc, ne, nm;
    char **chname = NULL;
    tnei_read_dataset("in/s019_r00418.txt", &val, &nt, &nc, &ne, &nm, &chname);
    assert(nm >= 2);
    char *capType = "R128"; /* !!! Should be a header field !!! */

    /* Define the imag size: */
    int32_t NX = 450;
    int32_t NY = 630;
    
    /* Read an image basis: */
    float_image_t **bas = tnei_read_basis_images("in/bdir/450_630_e128_gauss_i1_n1", ne, chname, NX, NY);
    
    fprintf(stderr, "--- testing {neuromat_eeg_image_make_idealized_scalp_mask} ---\n");
    r2_t ctr = (r2_t){{ (225.0/450.0)*NX, (315.0/630.0)*NY }};
    r2_t rad = (r2_t){{ (213.8/450.0)*NX, (256.5/630.0)*NY }};
    float_image_t *msk = neuromat_eeg_image_make_idealized_scalp_mask(NX, NY, &ctr, &rad);
    tnei_write_fni("msk", msk);
    
    fprintf(stderr, "--- testing {neuromat_eeg_image_compute_pot_field} ---\n");
    float_image_t *fld = float_image_new(1, NX, NY);
    int32_t it_x = (int32_t)(0.45*nt);  /* A random frame index. */
    neuromat_eeg_image_compute_pot_field(ne, val[it_x], bas, msk, fld);
    tnei_write_fni("fld", fld);
    
    fprintf(stderr, "--- testing {neuromat_eeg_image_electrodes_overlay} ---\n");
    float_image_t *cim = float_image_new(4, NX, NY);
    float_image_fill(cim, 0.000f);
    frgb_t fcfill = (frgb_t){{ 1.000f, 0.950f, 0.000f }};
    frgb_t fcdraw = (frgb_t){{ 0.000f, 0.350f, 1.000f }};
    double drad_el = 5.0;
    double hwd_el = 2.0;
    int32_t ie_spec = 3; /* Electrode to highlight. */
    int32_t ne_full = -1;
    r2_t *pos2D = NULL;
    char **chname_full = NULL;
    neuromat_eeg_geom_get_schematic_2D_points(capType, &ne_full, &chname_full, &pos2D);
    demand(ne_full == ne, "wrong electrode count");
    r2_t pos_img[ne]; /* Position of elecrodes on image. */
    neuromat_eeg_geom_map_many_disk_to_ellipse(ne, pos2D, &ctr, &rad, pos_img);
    float_image_t *elp = neuromat_eeg_image_electrodes_overlay(ne, pos_img, drad_el, hwd_el, ie_spec, &fcfill, &fcdraw, NX, NY);
    tnei_write_fni("elp", elp);
    
    fprintf(stderr, "--- testing {neuromat_eeg_image_paint_timeline_bars} ---\n");
    int32_t hwd_tb = 2;
    neuromat_eeg_marker_spec_t marker[nm];
    for (uint32_t im = 0;  im < nm; im++) 
      { marker[im].ic = ne + im;
        double ang = M_PI*((double)im)/(nm-1);
        marker[im].vref = 1.00;
        frgb_t fc = (frgb_t){{ (float)(0.5*(1 + cos(ang))), 1.000f, (float)(0.5*(1 - cos(ang))) }};
        double fcm = fmax(fc.c[0], fc.c[2]);
        fc.c[0] = (float)(fc.c[0]/fcm);
        fc.c[2] = (float)(fc.c[2]/fcm);
        marker[im].color = fc;
      }
    int32_t mkdots_xctr;
    double mkdots_rad;
    int32_t track_xlo, track_xsz;
    int32_t track_y[nm];
    int32_t slider_hh;
    float_image_t *trk = neuromat_eeg_image_make_time_tracks
      ( nt, nc, val, nm,  marker, hwd_tb, NX, &mkdots_xctr, &mkdots_rad, &track_xlo, &track_xsz, track_y, &slider_hh );
    int32_t slider_hw = 1;
    double tlo = 0.0;
    double thi = (double)nt;
    double t = 0.67*tlo + 0.33*thi;
    int32_t ylo = track_y[0];
    int32_t yhi = track_y[nm-1];
    fprintf(stderr, "track y = {%d..%d}\n", ylo, yhi);
    frgb_t fslid = (frgb_t){{ 1.000f, 0.000f, 0.500f }};
    neuromat_image_paint_slider(trk, slider_hw, slider_hh, tlo, t, thi, track_xlo, track_xsz, ylo, yhi, fslid);

    fprintf(stderr, "--- testing {neuromat_eeg_image_paint_marker_dots} ---\n");
    r2_t ctr_mark[nm];
    for (uint32_t im = 0;  im < nm; im++) 
      { ctr_mark[im] = (r2_t){{ (double)mkdots_xctr, (double)track_y[im] }}; }
    neuromat_eeg_image_paint_marker_dots(trk, nc, val[it_x], nm, marker, ctr_mark, mkdots_rad);
    tnei_write_fni("trk", trk);

    float_image_free(msk);
    for (uint32_t ie = 0;  ie < ne; ie++) { float_image_free(bas[ie]); }
    free(bas);
    for (uint32_t it = 0;  it < nt; it++) { free(val[it]); }
    free(val);
    float_image_free(cim);
    
    return 0;
  }
    
void tnei_write_fni(char *tag, float_image_t *img)
  { char *fname = NULL;
    char *fname = jsprintf("out/%s.fni", tag);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, img);
    fclose(wr);
    free(fname);
  }
 
float_image_t **tnei_read_basis_images(char *prefix, int32_t ne, char *chname[], int32_t NX, int32_t NY)
  {
    float_image_t **bas = notnull(malloc(ne*sizeof(float_image_t*)), "no mem");
    int32_t ie;
    for (ie = 0; ie < ne; ie++)
      { char *fname = NULL;
        char *fname = jsprintf("%s_%s.fni", prefix, chname[ie]);
        FILE *rd = open_read(fname, TRUE);
        float_image_t *bsi = float_image_read(rd);
        fclose(rd);
        demand(bsi->sz[0] == 1, "electrode position image must have a single channel");
        demand(bsi->sz[1] == NX, "electrode position image has wrong width");
        demand(bsi->sz[2] == NY, "electrode position image has wrong height");
        bas[ie] = bsi;
        free(fname);
      }
    return bas;
  }

void tnei_read_dataset(char *fname, double ***valP, int32_t *ntP, int32_t *ncP, int32_t *neP, int32_t *nmP, char ***chnameP)
  { 
    FILE *rd = open_read(fname, TRUE);
    
    /* Read the EEG header: */
    int32_t nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(rd, 20, 600, &nl);
    int32_t nc = h->nc;    /* Number of channels per data frame (including markers etc.). */
    int32_t ne = h->ne;    /* Number of electrodes. */
    int32_t nm = nc - ne;  /* Number of marker channel dots. */
    char **chname = h->chname;
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
    (*chnameP) = chname;
  }


#define PROG_NAME "test_neuromat_eeg_image"
#define PROG_DESC "Tests {neuromat_eeg_image.h}."
#define PROG_VERS "2021-08-24"

#define test_neuromat_eeg_image_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"
/* Last edited on 2021-08-25 03:18:37 by stolfi */

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
 
void tnei_read_dataset(char *fname, double ***valP, int32_t *ntP, int32_t *ncP, int32_t *neP, int32_t *nmP, char ***chnamesP);
float_image_t **tnei_read_basis_images(char *prefix, int32_t ne, char *chnames[], int32_t NX, int32_t NY);
void tnei_write_fni(char *tag, float_image_t *img);

int32_t main(int32_t argc, char **argv)
  {
    
    /* Read an EEG dataset: */
    double **val = NULL;
    int32_t nt, nc, ne, nm;
    char **chnames = NULL;
    tnei_read_dataset("in/s019_r00418.txt", &val, &nt, &nc, &ne, &nm, &chnames);
    assert(nm >= 2);

    /* Define the imag size: */
    int32_t NX = 450;
    int32_t NY = 630;
    
    /* Read an image basis: */
    float_image_t **bas = tnei_read_basis_images("in/bdir/450_630_e128_gauss_i1_n1", ne, chnames, NX, NY);
    
    fprintf(stderr, "--- testing {neuromat_eeg_image_make_idealized_scalp_mask} ---\n");
    r2_t ctr = (r2_t){{ 0.50*NX, 0.65*NY }};
    r2_t rad = (r2_t){{ 0.35*NX, 0.35*NY }};
    float_image_t *msk = neuromat_eeg_image_make_idealized_scalp_mask(NX, NY, &ctr, &rad);
    tnei_write_fni("msk", msk);
    
    fprintf(stderr, "--- testing {neuromat_eeg_image_compute_pot_field} ---\n");
    float_image_t *fld = float_image_new(1, NX, NY);
    int32_t it_x = (int32_t)(0.45*nt);  /* A random frame index. */
    neuromat_eeg_image_compute_pot_field(ne, val[it_x], bas, msk, fld);
    tnei_write_fni("fld", fld);
    
    fprintf(stderr, "--- testing {neuromat_eeg_image_paint_electrodes} ---\n");
    float_image_t *cim = float_image_new(4, NX, NY);
    float_image_fill(cim, 0.000f);
    frgb_t fcfill = (frgb_t){{ 1.000f, 0.950f, 0.000f }};
    frgb_t fcdraw = (frgb_t){{ 0.000f, 0.350f, 1.000f }};
    double drad_el = 5.0;
    double hwd_el = 2.0;
    int32_t ie_spec = 3; /* Electrode to highlight. */
    r2_t *pos2D = neuromat_eeg_geom_get_schematic_2D_points(ne); 
    r2_t pos_img[ne]; /* Position of elecrodes on image. */
    neuromat_eeg_geom_map_many_disk_to_ellipse(ne, pos2D, &ctr, &rad, pos_img);
    neuromat_eeg_image_paint_electrodes(ne, pos_img, drad_el, hwd_el, ie_spec, &fcfill, &fcdraw, cim);
    tnei_write_fni("elp", cim);
    
    fprintf(stderr, "--- testing {neuromat_eeg_image_paint_timeline_bars} ---\n");
    float_image_fill(cim, 0.000f);
    int32_t hwd_tb = 2;
    int32_t ic_mark[nm];
    frgb_t fmark_p[nm];
    frgb_t fmark_n[nm];
    for (int32_t im = 0; im < nm; im++) 
      { ic_mark[im] = ne + im;
        double ang = M_PI*((double)im)/(nm-1);
        fmark_p[im] = (frgb_t){{ (float)(0.5*(1 + cos(ang))), 1.000f, (float)(0.5*(1 - cos(ang))) }};
        fmark_n[im] = (frgb_t){{ (float)(0.5*(1 - cos(ang))), 0.000f, (float)(0.5*(1 + cos(ang))) }};
      }
    int32_t track_xlo, track_xhi;
    int32_t track_y[nm];
    neuromat_eeg_image_paint_timeline_bars(nt, nc, val, nm,  ic_mark, hwd_tb, fmark_p, cim, &track_xlo, &track_xhi, track_y);
    tnei_write_fni("trk", cim);

    fprintf(stderr, "--- testing {neuromat_eeg_image_paint_timeline_bars} ---\n");
    float_image_fill(cim, 0.000f);
    double rad_mark = 15.0;
    double vscale[nm];
    r2_t ctr_mark[nm];
    for (int32_t im = 0; im < nm; im++) 
      { ic_mark[im] = ne + im;
        vscale[im] = 1.00;
        ctr_mark[im] = (r2_t){{ 20, NY - 20 - (2*rad_mark + 10)*im }}; 
      }
    neuromat_eeg_image_paint_marker_dots(nc, val[it_x], nm, ic_mark, vscale, fmark_p, fmark_n, ctr_mark, rad_mark, cim);
    tnei_write_fni("mdt", cim);

    float_image_free(msk);
    for (int32_t ie = 0; ie < ne; ie++) { float_image_free(bas[ie]); }
    free(bas);
    for (int32_t it = 0; it < nt; it++) { free(val[it]); }
    free(val);
    float_image_free(cim);
    
    return 0;
  }
    
void tnei_write_fni(char *tag, float_image_t *img)
  { char *fname = NULL;
    asprintf(&fname, "out/%s.fni", tag);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, img);
    fclose(wr);
    free(fname);
  }
 
float_image_t **tnei_read_basis_images(char *prefix, int32_t ne, char *chnames[], int32_t NX, int32_t NY)
  {
    float_image_t **bas = notnull(malloc(ne*sizeof(float_image_t*)), "no mem");
    int32_t ie;
    for (ie = 0; ie < ne; ie++)
      { char *fname = NULL;
        asprintf(&fname, "%s_%s.fni", prefix, chnames[ie]);
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

void tnei_read_dataset(char *fname, double ***valP, int32_t *ntP, int32_t *ncP, int32_t *neP, int32_t *nmP, char ***chnamesP)
  { 
    FILE *rd = open_read(fname, TRUE);
    
    /* Read the EEG header: */
    int32_t nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(rd, 20, 600, &nl);
    int32_t nc = h->nc;    /* Number of channels per data frame (including markers etc.). */
    int32_t ne = h->ne;    /* Number of electrodes. */
    int32_t nm = nc - ne;  /* Number of marker channel dots. */
    char **chnames = h->chnames;
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
    (*chnamesP) = chnames;
  }


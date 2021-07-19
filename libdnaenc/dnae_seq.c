/* See dnae_seq.h */
/* Last edited on 2021-07-18 00:42:10 by jstolfi */

#define dnae_seq_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <fftw3.h>

#include <vec.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <jsfile.h>
#include <hermite3.h>
#include <conv_filter.h>

#include <msm_seq_desc.h>

#include <dnae_nucleic.h>
#include <dnae_sample.h>
#include <dnae_datum.h>

#include <dnae_seq.h>

/* IMPLEMENTATIONS */

dnae_seq_t dnae_seq_new(int n)
  { dnae_seq_t seq;
    seq.sd.id = dnae_seq_id_none;
    seq.sd.name = NULL;
    seq.sd.rev = FALSE;
    seq.sd.size = n;
    seq.sd.estep = 0;
    seq.sd.skip = 0;
    seq.cmt = NULL;
    int k;
    for (k = 0; k < dnae_CHANNELS; k++) { seq.sfac.f[k] = 1.0; }
    seq.dv = dnae_datum_vec_new(n);
    return seq;
  }

dnae_seq_t dnae_seq_from_datum_vec
  ( dnae_seq_id_t id, 
    char *name, 
    bool_t rev,
    int8_t estep, 
    int skip, 
    char *cmt, 
    dnae_datum_scale_t *sfac,
    dnae_datum_vec_t dv
  )
  { dnae_seq_t seq;
    seq.sd.id = id;
    seq.sd.name = name;
    seq.sd.rev = FALSE;
    seq.sd.size = dv.ne;
    seq.sd.estep = estep;
    seq.sd.skip = skip;
    seq.cmt = cmt;
    int k;
    for (k = 0; k < dnae_CHANNELS; k++) { seq.sfac.f[k] = sfac->f[k]; }
    seq.dv = dv;
    return seq;
  }
  
dnae_seq_t dnae_seq_from_nucleic_string(dnae_seq_id_t id, char *name, bool_t rev, char *cmt, char *bas)
  { demand(! rev, "reversion not supported yet");
    /* Convert nucleotide codes to numeric datums: */
    dnae_datum_vec_t dv = dnae_datum_vec_from_nucleic_string(bas); /* Should pass {rev}. */
    /* Compute standard deviation {rms} of signal for uniformly distributed bases: */
    double rms = dnae_NUCLEIC_RAW_SCALE;
    dnae_datum_scale_t sfac;
    int k; for (k = 0; k < dnae_CHANNELS; k++) { sfac.f[k] = rms; }
    /* Pack it all as a {dnae_seq_t}: */
    int8_t estep = 0;
    int skip = 0;
    dnae_seq_t seq = dnae_seq_from_datum_vec(id, name, rev, estep, skip, cmt, &sfac, dv); 
    return seq;
  }

dnae_seq_t dnae_seq_copy(dnae_seq_t *seq)
  { int nd = seq->sd.size;
    return dnae_seq_copy_sub(seq, 0, nd - 1); 
  }

dnae_seq_t dnae_seq_copy_sub(dnae_seq_t *seq, int ini, int fin)
  { /* Compute the new sequence's size and offset: */
    int size = (ini > fin ? 0 : fin - ini + 1);
    int skip = (size == 0 ? 0 : seq->sd.skip + ini);
    if (size > 0) { demand((0 <= ini) && (fin < seq->sd.size), "invalid index range"); }
    /* Create the descriptor and copy the datums: */
    dnae_seq_t scp;
    scp.sd.id = seq->sd.id;
    scp.sd.name = seq->sd.name;
    scp.sd.rev = seq->sd.rev;
    scp.sd.estep = seq->sd.estep;
    scp.sd.skip = skip;
    scp.sd.size = size;
    scp.cmt = seq->cmt;
    int i, k;
    for (k = 0; k < dnae_CHANNELS; k++) { scp.sfac.f[k] = seq->sfac.f[k]; }
    scp.dv = dnae_datum_vec_new(size);
    for (i = 0; i < size; i++) { scp.dv.e[i] = seq->dv.e[i + ini]; }
    return scp; 
  }

int dnae_seq_num_datums(dnae_seq_t *seq)
  { return seq->dv.ne; }
  
dnae_sample_enc_t *dnae_seq_get_sample_enc_address(dnae_seq_t *seq, int i, int k)
  { assert((k >= 0) && (k < dnae_CHANNELS));
    dnae_datum_t *dp = dnae_seq_get_datum_address(seq, i);
    return &(dp->c[k]);
  }
  
dnae_sample_enc_t dnae_seq_get_sample_enc(dnae_seq_t *seq, int i, int k)
  { dnae_sample_enc_t *sp = dnae_seq_get_sample_enc_address(seq, i, k);
    return *sp;
  }
  
void dnae_seq_set_sample_enc(dnae_seq_t *seq, int i, int k, dnae_sample_enc_t s)
  { dnae_sample_enc_t *sp = dnae_seq_get_sample_enc_address(seq, i, k);
    (*sp) = s;
  }

dnae_datum_t *dnae_seq_get_datum_address(dnae_seq_t *seq, int i)
  { int n = seq->dv.ne;
    
    assert((i >= 0) && (i < n));
    dnae_datum_t *dp = &(seq->dv.e[i]);
    return dp;
  }

dnae_datum_t dnae_seq_get_datum(dnae_seq_t *seq, int i)
  { dnae_datum_t *dp = dnae_seq_get_datum_address(seq, i);
    return *dp;
  }
  
void dnae_seq_set_datum(dnae_seq_t *seq, int i, dnae_datum_t d)
  { dnae_datum_t *dp = dnae_seq_get_datum_address(seq, i);
    (*dp) = d;
  }

dnae_datum_t dnae_seq_eval(dnae_seq_t *seq, double x, bool_t smooth)
  { 
    int nsmp = seq->dv.ne;
    int32_t size = seq->sd.size;
    assert(nsmp == size);
    
    /* Check range of position {x}: */
    demand((x >= 0) && (x <= size-1), "invalid sequence position");
    
    /* Linear interpolation for now. */
    /* !!! Change to a smoother interpolation !!! */
    int i0 = (int)floor(x); /* Index of preceding sample. */
    int i1 = (int)ceil(x);  /* Index of next sample. */
    assert((0 <= i0) && (i0 <= i1) && (i1 < nsmp));
    if (i0 == i1)
      { /* No interpolation needed. */
        return seq->dv.e[i0];
      }
    else 
      { /* Interpolation between datums at {i0} and {i1}: */
        double r = x - i0;
        assert((0.0 <= r) && (r <= 1.0));
        assert(nsmp >= 2);
        /* Get the two enclosing samples: */
        dnae_datum_t *f0 = &(seq->dv.e[i0]);
        dnae_datum_t *f1 = &(seq->dv.e[i1]);
        /* Get an extra two samples on each side if smooth interpolation: */
        dnae_datum_t *fmm = NULL, *fm = NULL, *fp = NULL, *fpp = NULL;
        if (smooth)
          { fmm = (i0 < 2 ? NULL : &(seq->dv.e[i0-2]));
            fm = (i0 < 1 ? NULL : &(seq->dv.e[i0-1]));
            fp = (i1 > nsmp-2 ? NULL : &(seq->dv.e[i1+1]));
            fpp = (i1 > nsmp-3 ? NULL : &(seq->dv.e[i1+2]));
          }
        
        /* Evaluate the interpolated datum {fx} channel by channel: */
        dnae_datum_t fx;
        int c;
        for (c = 0; c < dnae_CHANNELS; c++)
          { double sfac = seq->sfac.f[c]; /* Scale factor for channel {c}. */
            /* Get decoded values at {i0,i1}: */
            double v0 = dnae_sample_decode(f0->c[c], sfac);
            double v1 = dnae_sample_decode(f1->c[c], sfac);
            /* Interpolate: */
            double vx; /* Interpolated value. */
            if ((nsmp == 2) || (! smooth))
              { /* Linear interpolation:  */
                vx = (1-r)*v0 + r*v1;
              }
            else if (nsmp == 3)
              { /* Parabola:  */
                if (i0 == 0)
                  { assert(fp != NULL);
                    double vp = dnae_sample_decode(fp->c[c], sfac);
                    vx = 0.5*(r - 1)*(r - 2)*v0 - r*(r - 2)*v1 + 0.5*r*(r - 1)*vp;
                  }
                else
                  { assert(fm != NULL);
                    double vm = dnae_sample_decode(fp->c[c], sfac);
                    vx = 0.5*(r - 1)*(r - 2)*vm - r*(r - 2)*v0 + 0.5*r*(r - 1)*v1;
                  }
              }
            else
              { /* Cubic interpolation. */
                /* Extract the values {i0-2,i0-1,i1+1,i1+2} ix they exist and are needed: */
                double vmm = ( i1 == nsmp-1 ? dnae_sample_decode(fmm->c[c], sfac) : NAN );
                double vm = ( i0 > 0 ? dnae_sample_decode(fm->c[c], sfac) : NAN );
                double vp = ( i1 < nsmp-1 ? dnae_sample_decode(fp->c[c], sfac) : NAN );
                double vpp = ( i0 == 0 ? dnae_sample_decode(fpp->c[c], sfac) : NAN );
                
                /* Compute derivatives {d0,d1} at {i0,i1} from data points that exist: */
                double d0, d1;
                if (i0 == 0)
                  { assert(! isnan(vpp));
                    assert(! isnan(vp));
                    d0 = hermite3_estimate_deriv_0_1_3(v0, v1, vp, vpp);
                    d1 = hermite3_estimate_deriv_1_1_2(v0, v1, vp, vpp);
                  }
                else if (i1 == nsmp-1)
                  { assert(! isnan(vmm));
                    assert(! isnan(vm));
                    d0 = - hermite3_estimate_deriv_1_1_2(v1, v0, vm, vmm);
                    d1 = - hermite3_estimate_deriv_0_1_3(v1, v0, vm, vmm);
                  }
                else 
                  { assert(! isnan(vm));
                    assert(! isnan(vp));
                    if (i0 == 1)
                      { d0 = hermite3_estimate_deriv_1_1_2(vm, v0, v1, vp); }
                    else
                      { assert(! isnan(vmm));
                        d0 = hermite3_estimate_deriv_2_0_2(vmm, vm, v1, vp);
                      }
                    if (i1 == nsmp-2)
                      { d1 = - hermite3_estimate_deriv_1_1_2(vp, v1, v0, vm); }
                    else
                      { assert(! isnan(vpp));
                        d1 = hermite3_estimate_deriv_2_0_2(vm, v0, vp, vpp);
                      }
                  }
                /* Interpolate with cubic Hermite: */
                vx = hermite3_interp(v0, d0, v1, d1, 0.0, 1.0, r);
              }
            /* Re-encode and store in result: */
            fx.c[c] = dnae_sample_encode(vx, sfac);
          }
        return fx;
      }
  }

void dnae_seq_write(FILE *wr, dnae_seq_t *seq)
  { 
    /* Write the file header: */
    filefmt_write_header(wr, dnae_seq_type_name, dnae_seq_version);
    
    /* Write the comment, if any: */
    int ind = 0; /* Comment indentation. */
    if (seq->cmt != NULL) { filefmt_write_comment(wr, seq->cmt, ind, '|'); }
    
    /* Write the header fields: */
    fprintf(wr, "id = %d\n", seq->sd.id);
    fprintf(wr, "name = %s\n", seq->sd.name);
    fprintf(wr, "reversed = %c\n", ("FT"[seq->sd.rev]));
    fprintf(wr, "resampling =  %d\n", seq->sd.estep);
    fprintf(wr, "offset = %d\n", seq->sd.skip);
    fprintf(wr, "channels = %d\n", dnae_CHANNELS);
    fprintf(wr, "scale =");
    int k; for (k = 0; k < dnae_CHANNELS; k++)
      { fprintf(wr, " %24.16e", seq->sfac.f[k]); }
    fprintf(wr, "\n");
    fprintf(wr, "datums = %d\n", seq->dv.ne);
    
    /* Write the datum vector: */
    int ib;
    for (ib = 0; ib < seq->dv.ne; ib++)
      { dnae_datum_t *dp = &(seq->dv.e[ib]);
        int k;
        for (k = 0; k < dnae_CHANNELS; k++)
          { fprintf(wr, " %+6d", dp->c[k]); }
        fprintf(wr, "\n");
      }

    /* Write the file footer: */
    filefmt_write_footer(wr, dnae_seq_type_name);
    fflush(wr);
  }

void dnae_seq_write_named(char *fname, char *tag, char *ext, dnae_seq_t *seq)
  { FILE *wr = open_write_tag_ext(fname, tag, ext, TRUE);
    dnae_seq_write(wr, seq);
    if (wr != stdout) { fclose(wr); }
  }

dnae_seq_t dnae_seq_read(FILE *rd)
  { 
    /* Check and skip the file header: */
    filefmt_read_header(rd, dnae_seq_type_name, dnae_seq_version);
    
    /* Skip comment lines: */
    char *cmt = filefmt_read_comment(rd, '|');
    
    /* Read the header fields: */
    int id = nget_int(rd, "id"); fget_eol(rd);

    char *name = nget_string(rd, "name"); fget_eol(rd);

    bool_t rev = nget_bool(rd, "reversed"); fget_eol(rd);

    int estep = nget_int(rd, "resampling"); fget_eol(rd);
    demand(abs(estep) <= msm_seq_desc_estep_MAX, "bad resampling");

    int skip = nget_int(rd, "offset"); fget_eol(rd);

    int chns = nget_int(rd, "channels"); fget_eol(rd);
    demand(chns == dnae_CHANNELS, "wrong number of channels");
    
    dnae_datum_scale_t sfac;
    sfac.f[0] = nget_double(rd, "scale");
    int k; for (k = 1; k < dnae_CHANNELS; k++) { sfac.f[k] = fget_double(rd); }
    fget_eol(rd);
    
    int nsmp = nget_int(rd, "datums"); fget_eol(rd);
    demand(nsmp >= 0, "invalid datum count");
    
    /* Read the datum vector: */
    dnae_datum_vec_t dv = dnae_datum_vec_new(nsmp);
    int ib;
    for (ib = 0; ib < nsmp; ib++)
      { dnae_datum_t d;
        int k;
        for (k = 0; k < dnae_CHANNELS; k++)
          { int s = fget_int(rd);
            demand((s >= dnae_sample_enc_VALID_MIN) && (s <= dnae_sample_enc_VALID_MAX), "bad sample value"); 
            d.c[k] = (dnae_sample_enc_t)s;
          }
        fget_eol(rd);
        dv.e[ib] = d;
      }
    
    /* Check and skip file footer: */
    filefmt_read_footer(rd, dnae_seq_type_name);
    
    /* Create the {dnae_seq_t} record: */
    dnae_seq_t seq = dnae_seq_from_datum_vec(id, name, rev, (int8_t)estep, skip, cmt, &sfac, dv);
    return seq;
  }

dnae_seq_t dnae_seq_read_named(char *fname, char *tag, char *ext)
  { FILE *rd = open_read_tag_ext(fname, tag, ext, TRUE);
    dnae_seq_t s = dnae_seq_read(rd);
    if (rd != stdin) { fclose(rd); }
    return s;
  }

dnae_seq_t dnae_seq_read_from_nucleic_file(FILE *rd, dnae_seq_id_t id, char *name, bool_t rev)
  { /* Read letter sequence from file: */
    char *bas = NULL;
    char *cmt = NULL;
    dnae_nucleic_string_read(rd, &bas, &cmt);
    /* Pack it all up: */
    dnae_seq_t seq = dnae_seq_from_nucleic_string(id, name, rev, cmt, bas);
    free(bas);
    return seq;
  }
  
dnae_seq_t dnae_seq_read_from_nucleic_file_named(char *fname, char *tag, char *ext, dnae_seq_id_t id, char *name, bool_t rev)
  { /* Read letter sequence from file "{fname}{tag}{ext}": */
    char *bas = NULL;
    char *cmt = NULL;
    dnae_nucleic_string_read_named(fname, tag, ext, &bas, &cmt);
    /* Pack it all up: */
    dnae_seq_t seq = dnae_seq_from_nucleic_string(id, name, rev, cmt, bas);
    free(bas);
    return seq;
  }
  
void dnae_seq_free_datums(dnae_seq_t *seq)
  {
    free(seq->dv.e); 
  }
    
void dnae_seq_free(dnae_seq_t *seq)
  {
    dnae_seq_free_datums(seq);
    free(seq); 
  }
    
void dnae_seq_multi_free_datums(dnae_seq_t seq[], int maxLevel)
  { int level;
    for (level = 0; level <= maxLevel; level++)
      { dnae_seq_free_datums(&(seq[level])); }
  }

void dnae_seq_multi_free(dnae_seq_t *seq[], int maxLevel)
  { int level;
    for (level = 0; level <= maxLevel; level++)
      { dnae_seq_free(seq[level]); }
  }
   
dnae_seq_t dnae_seq_filter(dnae_seq_t *seq, double_vec_t *wtb, int8_t ek, char *wcmt)
  { /* Get filter kernel width {nw}: */
    int nw = wtb->ne;
    demand((nw % 2) == 1, "filter kernel width must be odd");
    int hw = (nw - 1)/2;  /* Filter kernel radius. */
    /* Decide parameters of filtered sequence (note: may be empty): */
    demand(ek >= 0, "{dnae_seq_filter} cannot interpolate");
    msm_seq_desc_t sdo = seq->sd; /* Old sequence descriptor. */
    msm_seq_desc_t sdn = msm_seq_desc_filter(&sdo, nw, ek); /* New sequence descriptor. */
    assert(sdn.estep - sdo.estep == ek);
    int32_t ndo = sdo.size;
    int32_t ndn = sdn.size;
    assert(ndo == seq->dv.ne);
    /* Allocate new datum vector {dvn}: */
    dnae_datum_vec_t dvn = dnae_datum_vec_new(ndn);
    dnae_datum_scale_t Sn;  /* Scale factors for encoding the new samples: */
    int k;
    if (ndn == 0)
      { /* Sequence {seq} is empty, output is empty. */
        /* Set {Sn} to something sensible, just in case: */
        double den = pow(2.0, 0.5*ek);
        for (k = 0; k < dnae_CHANNELS; k++) { Sn.f[k] = seq->sfac.f[k]/den; }
        /* Vector {dvn} is empty. */
      }
    else
      { /* Some samples of {sdo} are aligned with some samples of {sdn}. */
        /* The assertions marked "OC" are output conditions of {msm_seq_desc_filter}. */
        int stpo = (1 << ek);  /* Step in {sdo} between aligned samples. */
        int32_t iskpo = sdn.skip*stpo - sdo.skip; /* Samples skipped at start of old sequence. */
        int32_t fskpo = ndo - iskpo - 1 - (ndn - 1)*stpo; /* Samples skipped at end of old sequence. */
        /* Old samples lost must be enough to fit the filtering kernel: */
        assert(iskpo >= hw); /* OC */
        assert(fskpo >= hw); /* OC */
        /* There should be some old samples left: */
        assert(iskpo + fskpo < ndo); /* OC */

        /* Process datum vector, one channel at a time; */
        double_vec_t vo = double_vec_new(ndo); /* Decoded samples of old sequence. */
        double_vec_t vn = double_vec_new(ndn); /* Decoded samples of new sequence. */
        for (k = 0; k < dnae_CHANNELS; k++)
          { /* Extract and decode the old samples {seq.dv}: */
            double Sok = seq->sfac.f[k];
            int i;
            for (i = 0; i < ndo; i++)
              { dnae_sample_enc_t s = dnae_seq_get_sample_enc(seq, i, k);
                vo.e[i] = dnae_sample_decode(s, Sok); 
              }

            /* Filter the old samples to produce the new samples {vn}, without interpolation: */
            conv_filter(ndo, vo.e, iskpo, stpo, nw, wtb->e, ndn, vn.e);
            
            /* Compute the new scale factor {Sn.f[k]}: */
            double sum_sq = 0;
            for (i = 0; i < ndn; i++) { double v = vn.e[i]; sum_sq += v*v; }
            Sn.f[k] = sqrt(sum_sq/ndn);
            /* fprintf(stderr, "  channel %d sum_sq = %14.8f  sfac = %12.8f\n", k, sum_sq, Sn.f[k]); */

            /* Encode the new samples and store them into the new datum vector: */
            for (i = 0; i < ndn; i++)
              { dnae_sample_enc_t s = dnae_sample_encode(vn.e[i], Sn.f[k]);
                dnae_datum_t *dp = &(dvn.e[i]);
                dp->c[k] = s; 
              }
          }
        free(vo.e);
        free(vn.e);
      }
      
    /* Append the weight table comment to the sequence's comment string: */
    if (wcmt == NULL) { wcmt = "unnamed table"; }
    char *old_cmt = (seq->cmt != NULL ? seq->cmt : "");
    char *sep = (strlen(old_cmt) == 0 ? "" : " ");
    char *new_cmt = NULL;
    asprintf(&(new_cmt), "%s%sfiltered with %s (window width %d)\n", old_cmt, sep, wcmt, nw);

    /* Pack the filtered datum vector as a {dnae_seq_t}: */
    dnae_seq_t new = dnae_seq_from_datum_vec
      ( sdn.id, sdn.name, sdn.rev, sdn.estep, sdn.skip,  
        new_cmt, &Sn, dvn
      );
    return new;
  }

dnae_seq_t dnae_seq_interpolate(dnae_seq_t *seq, int8_t ek)
  { /* Decide parameters of filtered sequence (note: may be empty): */
    demand(ek <= 0, "{dnae_seq_interpolate} cannot downsample");
    msm_seq_desc_t sdo = seq->sd; /* Old sequence descriptor. */
    msm_seq_desc_t sdn = msm_seq_desc_resample(&sdo, ek); /* New sequence descriptor. */
    assert(sdn.estep - sdo.estep == ek);
    int32_t ndo = sdo.size;
    int32_t ndn = sdn.size;
    assert(ndo == seq->dv.ne);
    int stpn = (1 << (-ek));  /* Step in {sdn} between aligned samples. */
    /* Allocate new datum vector {dvn}: */
    dnae_datum_vec_t dvn = dnae_datum_vec_new(ndn);
    dnae_datum_scale_t Sn;  /* Scale factors for encoding the new samples: */
    int k;
    if (ndn == 0)
      { /* Sequence {seq} is empty, output is empty. */
        /* Set {Sn} to something sensible, just in case: */
        double den = pow(2.0, 0.5*ek);
        for (k = 0; k < dnae_CHANNELS; k++) { Sn.f[k] = seq->sfac.f[k]/den; }
        /* Vector {dvn} is empty. */
      }
    else
      { /* Some samples of {sdn} are aligned with some samples of {sdo}. */
        /* The assertions marked "OC" are output conditions of {msm-seq_desc_interpolate}. */
        if (stpn > 1)
          { /* Process datum vector, one channel at a time; */
            double_vec_t vo = double_vec_new(ndo); /* Decoded samples of old sequence. */
            double_vec_t vn = double_vec_new(ndn); /* Decoded samples of new sequence. */
            for (k = 0; k < dnae_CHANNELS; k++)
              { /* Extract and decode the old samples {seq.dv}: */
                double Sok = seq->sfac.f[k];
                int i;
                for (i = 0; i < ndo; i++)
                  { dnae_sample_enc_t s = dnae_seq_get_sample_enc(seq, i, k);
                    vo.e[i] = dnae_sample_decode(s, Sok); 
                  }

                /* Interpolate the old samples to produce the new samples {vn}: */
                hermite3_subsample(ndo, vo.e, NULL, stpn, ndn, vn.e);

                /* Compute the new scale factor {Sn.f[k]}: */
                double sum_sq = 0;
                for (i = 0; i < ndn; i++) { double v = vn.e[i]; sum_sq += v*v; }
                Sn.f[k] = sqrt(sum_sq/ndn);
                /* fprintf(stderr, "  channel %d sum_sq = %14.8f  sfac = %12.8f\n", k, sum_sq, Sn.f[k]); */

                /* Encode the new samples and store them into the new datum vector: */
                for (i = 0; i < ndn; i++)
                  { dnae_sample_enc_t s = dnae_sample_encode(vn.e[i], Sn.f[k]);
                    dnae_datum_t *dp = &(dvn.e[i]);
                    dp->c[k] = s; 
                  }
              }
            free(vo.e);
            free(vn.e);
          }
      }
      
    /* Append the interpolation comment to the sequence's comment string: */
    char *old_cmt = (seq->cmt != NULL ? seq->cmt : "");
    char *sep = (strlen(old_cmt) == 0 ? "" : " ");
    char *new_cmt = NULL;
    asprintf(&(new_cmt), "%s%sinterpolated with step {1/%d}\n", old_cmt, sep, stpn);

    /* Pack the filtered datum vector as a {dnae_seq_t}: */
    dnae_seq_t new = dnae_seq_from_datum_vec
      ( sdn.id, sdn.name, sdn.rev, sdn.estep, sdn.skip,  
        new_cmt, &Sn, dvn
      );
    return new;
    
  }

void dnae_seq_multi_filter
  ( dnae_seq_t *s, 
    int maxLevel, 
    double_vec_t *wtb0,
    char *wname0,
    double_vec_t *wtb1,
    char *wname1,
    int8_t ek0,
    dnae_seq_t sr[]
  )
  { int k; /* Level of hierarchy. */
    for (k = 0; k <= maxLevel; k++)
      { if (k == 0)
          { sr[k] = dnae_seq_copy(s); }
        else
          { double_vec_t *wtb = (k == 1 ? wtb0 : wtb1);
            char *wname = (k == 1 ? wname0 : wname1);
            int8_t ekk = (int8_t)(k == 1 ? ek0 : 1);    /* Exponent of step for upsampling/downsampling. */
            int8_t ekf = (int8_t)(ekk >= 0 ? ekk : 0);  /* Exponent of step for downsampling. */
            sr[k] = dnae_seq_filter(&(sr[k-1]), wtb, ekf, wname);
            if (ekk < 0) 
              { dnae_seq_t tmp = sr[k];
                sr[k] = dnae_seq_interpolate(&tmp, ekk);
                dnae_seq_free_datums(&tmp);
              }
          }
      }
  }
 

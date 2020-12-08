/* See msm_seq_desc.h */
/* Last edited on 2018-03-04 22:58:21 by stolfilocal */

#define msm_seq_desc_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include <vec.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <msm_basic.h>

#include <msm_seq_desc.h>

msm_seq_desc_t msm_seq_desc_make
  ( msm_seq_id_t id,  /* Internal identifier of the sequence. */
    char *name,       /* External identifier of the sequence. */
    bool_t rev,       /* TRUE if the original sequence was reversed. */
    int32_t size,     /* Number of distinct values in the sequence. */
    int8_t estep,     /* Each sampling step is {2^estep} samples of original; may be negative. */
    int32_t skip      /* First sample matches sample {skip*(2^estep)} of original seq. */
  )
  { demand(size >= 0, "invalid {size}");
    demand(abs((int)estep) <= msm_seq_desc_estep_MAX, "invalid {estep}");
    msm_seq_desc_t s = (msm_seq_desc_t)
      { .id = id,
        .name = name,
        .rev = rev,
        .size = size,
        .estep = estep,
        .skip = skip
      };
    return s;
  }

bool_t msm_seq_desc_same_orig_seq(msm_seq_desc_t *sa, msm_seq_desc_t *sb, bool_t die)
  { if (sa->id != sb->id)
      { fail_test(die, "sequence {id}s differ"); }
    if ((sa->name != sb->name) && (strcmp(sa->name,sb->name) != 0))
      { fail_test(die, "sequence {name}s differ"); }
    if (sa->rev != sb->rev)
      { fail_test(die, "sequence {rev}s differ"); }
    return TRUE;
  }

bool_t msm_seq_desc_equal(msm_seq_desc_t *sa, msm_seq_desc_t *sb, bool_t die)
  { 
    if (sa->id != sb->id)
      { fail_test(die, "sequence {ids differ"); }
    if ((sa->name != sb->name) && (strcmp(sa->name,sb->name) != 0))
      { fail_test(die, "sequence {names differ"); }
    if (sa->rev != sb->rev)
      { fail_test(die, "sequence {rev}s differ"); }
    if (sa->size != sb->size)
      { fail_test(die, "sequence {size}s differ"); }
    if (sa->size > 0)
      { if (sa->estep != sb->estep)
          { fail_test(die, "sequence {estep} fields differ"); }
        if (sa->skip != sb->skip)
          { fail_test(die, "sequence {skip} fields differ"); }
      }
    return TRUE;
  }

void msm_seq_desc_adjust_skip_and_size(int64_t *skipP, int64_t *sizeP, int64_t m);
  /* Rounds {*skipP} up and {*sizeP} down so that {*skipP} and {(*sizeP)-1} 
    are multiples of {m}.  Requires {*sizeP} to be positive on
    input.  However, if the adjustment would result in {*sizeP} negative or zero,
    sets both to zero. */
      
void msm_seq_desc_adjust_skip_and_size(int64_t *skipP, int64_t *sizeP, int64_t m)
  {
    int64_t skip = (*skipP);
    int64_t size = (*sizeP);
    assert(size > 0);
    int64_t xitrim = iceil(skip, m) - (int64_t)skip;
    if (xitrim != 0)
      { /* Adjust {skip,size} for first sample: */
        assert(xitrim > 0);
        size = (size <= xitrim ? 0 : size - xitrim);
        skip = (size <= 0 ? 0 : skip + xitrim); 
      }
    int64_t steps = (size == 0 ? 0 : size - 1); /* Number of steps. */
    int64_t xftrim = (steps % m);
    if (xftrim != 0)
      { /* Adjust {size} for last sample: */
        assert(xftrim > 0);
        assert(xftrim < size); /* If {steps} is not 0 mod {m}, at least 1 sample will remain. */ 
        size -= xftrim; 
      }
    (*skipP) = skip;
    (*sizeP) = size;
  }

msm_seq_desc_t msm_seq_desc_trim(msm_seq_desc_t *s, int32_t itrim, int32_t ftrim)
  {
    /* Expand {s.size} and {s.skip} to 64 bits: */
    int64_t ssize = (int64_t)s->size;
    int64_t sskip = (int64_t)s->skip;
    int64_t ntrim = (int64_t)itrim + (int64_t)ftrim; /* Total samples trimmed. */
    /* Compute the parameters of the {t} sequence in 64 bits to avoid overflow: */
    int64_t tsize = (ssize <= ntrim ? 0 : ssize - ntrim);
    int64_t tskip = (tsize <= 0 ? 0 : sskip + (int64_t)itrim);
    /* Check for overflow: */
    demand((int64_t)((int32_t)tsize) == tsize, "overflow in {size}");
    demand((int64_t)((int32_t)tskip) == tskip, "overflow in {skip}");
    /* Pack and ship: */
    msm_seq_desc_t t = (msm_seq_desc_t) 
      { .id = s->id,
        .name = s->name,
        .rev = s->rev,
        .size = (int32_t)tsize,
        .estep = s->estep,
        .skip = (int32_t)tskip
      };
    return t;
  }

msm_seq_desc_t msm_seq_desc_resample(msm_seq_desc_t *s, int8_t ek)
  { 
    msm_seq_desc_t t = (*s);
    
    /* Trivial case: */
    if (ek == 0) { return t; }
    
    /* Compute the {estep} of {t}: */
    int testep = (int)s->estep + (int)ek;
    demand(abs(testep) <= msm_seq_desc_estep_MAX, "overflow in {estep}");
    t.estep = (int8_t)testep;
    
    /* Expand {s.size} and {s.skip} to 64 bits to avoid overflow: */
    int64_t tskip = (int64_t)s->skip;
    int64_t tsize = (int64_t)s->size;
    if (tsize > 0)
      { if (ek > 0)
          { /* Downsampling by {2^ek}: */
            /* Trim {s} implicitly as needed, adjusting {tsize} and {tskip}: */
            int64_t m = (1 << ek);
            msm_seq_desc_adjust_skip_and_size(&tskip, &tsize, m);
            /* Divide {tsize} and {tskip-1} by {2^ek}: */
            if (tsize > 0)
              { tskip = tskip/m;
                tsize = (tsize - 1)/m + 1;
              }
          }
        else
          { /* Upsampling by {2^{-ek}}, no trimming needed: */
            int64_t m = (1 << (-ek));
            /* fprintf(stderr, "   upsampling by %d\n", m); */
            tsize = (tsize - 1)*m + 1;
            tskip = tskip*m;
          }
      }
      
    /* Check for 32-bit overflow: */
    demand((int64_t)((int32_t)tsize) == tsize, "overflow in {size}");
    demand((int64_t)((int32_t)tskip) == tskip, "overflow in {skip}");
    /* Pack and ship: */
    t.size = (int32_t)tsize;
    t.skip = (int32_t)tskip;
    return t;
  }

msm_seq_desc_t msm_seq_desc_filter(msm_seq_desc_t *sd, int nw, int8_t ek)
  { demand((nw > 0) && (nw % 2 == 1), "main filter width must be odd");
    int trim = (nw - 1)/2;
    msm_seq_desc_t tsd = msm_seq_desc_trim(sd, trim, trim);
    msm_seq_desc_t rsd = msm_seq_desc_resample(&tsd, ek);
    return rsd;
  }

double msm_seq_desc_map_index_to_orig_seq(double ix, msm_seq_desc_t *s)
  {
    double step = pow(2.0, s->estep);
    return (((double)ix) + ((double)s->skip)) * step;
  }

double msm_seq_desc_map_index_from_orig_seq(double ix0, msm_seq_desc_t *s)
  { 
    double step = pow(2.0, s->estep);
    return ix0/step - ((double)s->skip);
  }
 
double msm_seq_desc_map_index(double ixa, msm_seq_desc_t *sa, msm_seq_desc_t *sb)
  {
    /* Must be the same original sequence: */
    (void)msm_seq_desc_same_orig_seq(sa, sb, TRUE);
    double ix0 = msm_seq_desc_map_index_to_orig_seq(ixa, sa);
    double ixb = msm_seq_desc_map_index_from_orig_seq(ix0, sb);
    return ixb;
  }

void msm_seq_desc_throw_segment(int32_t size, int32_t n, int32_t *iniP, int32_t *finP)
  { demand (n <= size, "sequence is too small");
    if (n == size) 
      { (*iniP) = 0;
        (*finP) = size - 1;
      }
    else
      { int32_t ix = int32_abrandom(0, size - n); 
        int32_t fx = ix + n - 1;
        assert((ix >= 0) && (ix < size));
        assert(ix <= fx); 
        assert(fx < size);
        (*iniP) = ix;
        (*finP) = fx;
      }
  }

void msm_seq_desc_write
  ( FILE *wr, 
    char *pre, 
    msm_seq_desc_t *s, 
    int idSize, 
    int nameSize, 
    int indexSize,
    char *suf
  )
  { if ((pre != NULL) & ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%*d", idSize, s->id);
    fprintf(wr, " %-*s", nameSize, s->name);
    fprintf(wr, " %c", (s->rev ? 'T' : 'F'));
    fprintf(wr, " %*d", indexSize, s->size);
    fprintf(wr, (s->estep == 0 ? " %3d" : " %+3d"), s->estep);
    fprintf(wr, " %*d", indexSize, s->skip);
    if ((suf != NULL) & ((*suf) != 0)) { fputs(suf, wr); }
  }

msm_seq_desc_t msm_seq_desc_read(FILE *rd, char *pre, char *suf)
  { if ((pre != NULL) & ((*pre) != 0))
      { fget_skip_spaces(rd); fget_match(rd, pre); }
    msm_seq_desc_t s = (msm_seq_desc_t)
      { .id = (msm_seq_id_t)fget_int(rd),
        .name = fget_string(rd),
        .rev = fget_bool(rd),
        .size = (int32_t)fget_int(rd),
        .estep = (int8_t)fget_int(rd),
        .skip = (int32_t)fget_int(rd)
      };
    (void)msm_seq_desc_is_valid(&s, TRUE);
    if ((suf != NULL) & ((*suf) != 0)) 
      { fget_skip_spaces(rd); fget_match(rd, suf); }
    return s;
  }

bool_t msm_seq_desc_is_valid(msm_seq_desc_t *s, bool_t die)
  { 
    if (s->id < 0)            
      { fail_test(die, "invalid {id}"); }
    if (s->name == NULL)      
      { fail_test(die, "null {name}"); }
    if (strlen(s->name) == 0) 
      { fail_test(die, "empty {name}"); }
    if (s->size < 0)          
      { fail_test(die, "invalid {size}"); }
    if (abs((int)s->estep) > msm_seq_desc_estep_MAX) 
      { fail_test(die, "invalid {estep}"); }
    /* Note that {s->size} may be negative. */
    if ((s->size > 0) && (s->estep < 0))
      { int32_t m = (1 << (- s->estep));
        if ((s->skip % m) != 0) 
          { fail_test(die, "{skip} is not multiple of {2^{-estep}}"); }
        if (((s->size - 1) % m) != 0) 
          { fail_test(die, "{size-1} is not multiple of {2^{-estep}}"); }
      }
    return TRUE;
  }

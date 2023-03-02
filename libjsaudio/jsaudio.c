/* See jsaudio.h */
/* Last edited on 2023-03-02 12:30:34 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsaudio.h>

/* IMPLEMENTATIONS */

sound_t jsa_allocate_sound(int32_t nc, int32_t ns)
  { sound_t s;
    demand (nc >= 1, "invalid number of channels");
    s.nc = nc;
    s.ns = ns;
    s.sv = malloc(nc*sizeof(float *));
    assert(s.sv != NULL);
    for (int32_t ic = 0; ic < nc; ic++)
      { s.sv[ic] = malloc(ns*sizeof(double));
        assert(s.sv[ic] != NULL);
        for (int32_t i = 0; i < ns; i++) { s.sv[ic][i] = 0.0; }
      }
    return s;
  }

sound_t jsa_copy_sound(sound_t *s, int32_t ini, int32_t ns)
  { assert(ns >= 0);
    assert((ini >= 0) && (ini < s->ns)); 
    sound_t r = jsa_allocate_sound(s->nc, ns);
    r.fsmp = s->fsmp;
    r.nc = s->nc;
    r.ns = ns;
    for (int32_t ic = 0; ic < s->nc; ic++)
      { for (int32_t i = 0; i < ns; i++) 
         { r.sv[ic][i] = s->sv[ic][i + ini]; }
      }
    return r;
  }

void jsa_add_sound(sound_t *s, int32_t sskip, sound_t *r, int32_t rskip, int32_t ns)
  { int32_t nc = (s->nc < r->nc ? s->nc : r->nc); 
    for (int32_t ic = 0; ic < nc; ic++)
      { for (int32_t i = 0; i < ns; i++) 
          { int32_t si = i + sskip;
            int32_t ri = i + rskip;
            if ((si >= 0) && (si < s->ns) && (ri >= 0) && (ri < r->ns))
              { s->sv[ic][ri] += r->sv[ic][si]; }
          }
       }
  }

void jsa_free_sound(sound_t *s)
  { if (s->sv != NULL) 
      { for (int32_t ic = 0; ic < s->nc; ic++)
          { if (s->sv[ic] != NULL) { free(s->sv[ic]); s->sv[ic] = NULL; } }
        free(s->sv); s->sv = NULL; 
      }
  }


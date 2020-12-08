/* See jsaudio.h */
/* Last edited on 2006-10-29 00:36:42 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <jsaudio.h>

/* IMPLEMENTATIONS */

sound_t jsa_allocate_sound(int nc, int ns)
  { sound_t s;
    s.sv = malloc(nc*sizeof(float *));
    assert(s.sv != NULL);
    s.nc = nc;
    int c;
    for (c = 0; c < nc; c++)
      { s.sv[c] = malloc(ns*sizeof(double));
        assert(s.sv[c] != NULL);
        int i;
        for (i = 0; i < ns; i++) { s.sv[c][i] = 0.0; }
      }
    return s;
  }

sound_t jsa_copy_sound(sound_t *s, int ini, int ns)
  { assert(ns >= 0);
    assert((ini >= 0) && (ini < s->ns)); 
    sound_t r = jsa_allocate_sound(s->nc, ns);
    r.fsmp = s->fsmp;
    r.nc = s->nc;
    r.ns = ns;
    int i;
    for (i = 0; i < ns; i++) 
      { int c;
        for (c = 0; c < s->nc; c++)
          { r.sv[c][i] = s->sv[c][i + ini]; }
      }
    return r;
  }

void jsa_add_sound(sound_t *s, int sskip, sound_t *r, int rskip, int ns)
  { int nc = (s->nc < r->nc ? s->nc : r->nc); 
    int i;
    for (i = 0; i < ns; i++) 
      { int si = i + sskip;
        int ri = i + rskip;
        if ((si >= 0) && (si < s->ns) && (ri >= 0) && (ri < r->ns))
          { int c;
            for (c = 0; c < nc; c++)
              { s->sv[c][ri] += r->sv[c][si]; }
          }
      }
  }

void jsa_free_sound(sound_t *s)
  { if (s->sv != NULL) 
      { int c;
        for (c = 0; c < s->nc; c++)
          { if (s->sv[c] != NULL) { free(s->sv[c]); s->sv[c] = NULL; } }
        free(s->sv); s->sv = NULL; 
      }
  }


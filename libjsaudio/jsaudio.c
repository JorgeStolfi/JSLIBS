/* See jsaudio.h */
/* Last edited on 2024-12-21 03:21:40 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsaudio.h>

/* IMPLEMENTATIONS */

jsaudio_t jsaudio_allocate_sound(uint32_t nc, uint32_t ns)
  { jsaudio_t s;
    demand (nc >= 1, "invalid number of channels");
    s.nc = nc;
    s.ns = ns;
    s.sv = talloc(nc, double*);
    assert(s.sv != NULL);
    for (int32_t ic = 0;  ic < nc; ic++)
      { s.sv[ic] = talloc(ns, double);
        assert(s.sv[ic] != NULL);
        for (int32_t i = 0;  i < ns; i++) { s.sv[ic][i] = 0.0; }
      }
    return s;
  }

jsaudio_t jsaudio_copy_sound(jsaudio_t *s, uint32_t ini, uint32_t ns)
  { assert(ini < s->ns); 
    jsaudio_t r = jsaudio_allocate_sound(s->nc, ns);
    r.fsmp = s->fsmp;
    r.nc = s->nc;
    r.ns = ns;
    for (int32_t ic = 0;  ic < s->nc; ic++)
      { for (int32_t i = 0;  i < ns; i++) 
         { r.sv[ic][i] = s->sv[ic][i + (int32_t)ini]; }
      }
    return r;
  }

void jsaudio_add_sound(jsaudio_t *s, uint32_t sskip, jsaudio_t *r, uint32_t rskip, uint32_t ns)
  { uint32_t nc = (s->nc < r->nc ? s->nc : r->nc); 
    for (int32_t ic = 0;  ic < nc; ic++)
      { for (int32_t i = 0;  i < ns; i++) 
          { uint32_t si = (uint32_t)i + sskip;
            uint32_t ri = (uint32_t)i + rskip;
            if ((si < s->ns) && (ri < r->ns))
              { s->sv[ic][ri] += r->sv[ic][si]; }
          }
       }
  }

void jsaudio_free_sound(jsaudio_t *s)
  { if (s->sv != NULL) 
      { for (int32_t ic = 0;  ic < s->nc; ic++)
          { if (s->sv[ic] != NULL) { free(s->sv[ic]); s->sv[ic] = NULL; } }
        free(s->sv); s->sv = NULL; 
      }
  }


/* See {float_image_trim.h}. */
/* Last edited on 2008-05-25 03:23:51 by stolfi */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
 
#include <bool.h>
#include <affirm.h>
#include <irange.h>
#include <frgb.h>
#include <float_image.h>

#include <float_image_trim.h>

/* INTERNAL PROTOTYPES */

void ftrm_exclude_edge(irange_t OBox[], int e);
  /* Shrinks the rectangle {OBox} by removing one row or column of pixels
     along edge number {e}. A no-op if {OBox} is empty. */

void ftrim_get_edge(irange_t OBox[], int e, irange_t *MBox);
  /* Sets {*MBox} to the rectangle containing the single
    row or column of pixels that lies along edge {e}
    of the rectangle {OBox}. */ 

bool_t ftrim_empty_box(irange_t OBox[]);
  /* TRUE iff the rectangle {OBox} is empty (has zero width or height). */

/* IMPLEMENTATIONS */

void float_image_trim_background(float_image_t *I, double noise, double except, irange_t OBox[], bool_t paddable[], frgb_t padcolor[])
  {
    /* Get image dimensions: */
    int NC, IX, IY;
    float_image_get_size(I, &NC, &IX, &IY);
    
    /* Start with the whole image: */
    OBox[0] = (irange_t){{0, IX}};
    OBox[1] = (irange_t){{0, IY}};
    
    /* Determine whether the image can be padded along each edge: */
    int e;
    for (e = 0; e < 4; e++)
      { irange_t MBox[2]; /* Row or column of pixels along edge {e}. */
        ftrim_get_edge(OBox, e, MBox);
        padcolor[e] = float_image_dominant_color(I, MBox, noise);
        paddable[e] = FALSE; /* For now, may be set TRUE later. */
      }
      
    /* Trim paddable edges: */
    bool_t try = TRUE;  /* FALSE when margins have been completely trimmed. */
    while ((! ftrim_empty_box(OBox)) && try)
      { try = FALSE; /* Unless we succeed: */
        irange_t MBox[2]; /* First and last pixels of edge {e}. */
        ftrim_get_edge(OBox, e, MBox);
        if (float_image_is_uniform(I, MBox, &(padcolor[e]), noise, except))
          { /* Edge {e} of {O} is part of margin: */
            paddable[e] = TRUE; 
            /* Remove one row/column of pixels along that edge: */
            ftrm_exclude_edge(OBox, e);
            /* We haven't finished yet. */
            try = TRUE; 
          }
      }
      
    /* Invalidate {padcolor[e]} for all edges {e} that were not trimmed: */
    for (e = 0; e < 4; e++)
      { if (! paddable[e]) { padcolor[e] = (frgb_t){{-1,-1,-1}}; } }
  }
  

/* Last edited on 2014-09-01 23:40:18 by stolfilocal */
/* See {dna_seq_view_tools.h} */

#define _GNU_SOURCE
#include <assert.h>

#include <GL/glu.h>
#include <GL/glut.h>

#include <jsrandom.h>
#include <frgb_ops.h>

#include <dna_seq_view_tools.h>

void dna_seq_view_tools_draw_tetrahedron(frgb_t *color)
  {
    /* Vertices: */
    r3_t v[4] = {
      (r3_t) {{+1,+1,+1}},
      (r3_t) {{+1,-1,-1}},
      (r3_t) {{-1,+1,-1}},
      (r3_t) {{-1,-1,+1}}
    } ;
    /* Set surface finish: */
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(color->c[0],color->c[1],color->c[2]);
    /* Draw the edges: */
    glBegin(GL_LINES);
    int i,j;
    for(i = 0; i < 3; i++)
      { GLfloat xi = (GLfloat)(v[i].c[0]);
        GLfloat yi = (GLfloat)(v[i].c[1]);
        GLfloat zi = (GLfloat)(v[i].c[2]);
        for(j = i+1; j < 4; j++)
          { GLfloat xj = (GLfloat)(v[j].c[0]);
            GLfloat yj = (GLfloat)(v[j].c[1]);
            GLfloat zj = (GLfloat)(v[j].c[2]);
            glVertex3f(xi, yi, zi);
            glVertex3f(xj, yj, zj);
          }
      }
    glEnd();
  }

void dna_seq_view_tools_draw_sequence
  ( dnae_seq_t *z, 
    double magnify, 
    r3_t *pert, 
    double radius, 
    int nst,
    frgb_t *color, 
    int ini, 
    int fin
  )
  {
    int nsmp = dnae_seq_num_datums(z);
    assert(dnae_CHANNELS == 3);
    assert(nst >= 1);
   
    frgb_t effcolor = (frgb_t){{ 0.8f, 0.8f, 0.8f }}; /* Color of effaced parts of sequence. */
    
    /* The output is the union of segments and original datum balls.
      Each segment consists one or more consecutive interpolated datum
      balls with the same texture, and some connecting rods.  A rod
      is visible iff both ends are visible.  Rods that end at original
      datums are trimmed.
    */

    int ne = -1;    /* Number of items (rods, balls) in the current segment; {-1} between segments. */
    bool_t seg_vis; /* If {ne>0}: true if current segment is fully visible, false if effaced. */
    
    auto void start_segment(void);
      /* Start a segment of a chaing with uniform visibility. */
    
    auto void extend_segment(bool_t vis);
      /* To be called before appending another item to the current segment. */
    
    auto void end_segment(void);
      /* End previous segment, which had visibility {vis}. */
    
    void start_segment(void)
      { assert(ne == -1);
        ne = 0;
        /* {seg_vis} is indeterminate. */
      }
    
    void extend_segment(bool_t vis)
      { assert(ne >= 0);
        if (ne == 0)
          { 
            /* Set surface finish: */
            if (vis)
              { /* Opaque with given color: */
                glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
                glColor3f(color->c[0],color->c[1],color->c[2]); 
                glEnable(GL_COLOR_MATERIAL);
              }
            else
              { /* Should set transparency: */
                glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
                glColor3f(effcolor.c[0],effcolor.c[1],effcolor.c[2]); 
                glEnable(GL_COLOR_MATERIAL);
              }
            seg_vis = vis;
          }
        else
          { /* Items in segment should have consistent visibility: */
            assert(vis == seg_vis);
          }
        ne++;
      }
    
    void end_segment(void)
      { assert(ne >= 0);
        if (ne > 0)
          { glDisable(GL_COLOR_MATERIAL); }
        ne = -1;
      }
    
    GLUquadricObj* quad =  gluNewQuadric();
    double datrad = radius; /* Radius of datum balls. */
    double cylrad = 0.30*datrad; /* Radius of connecting cylinders. */

    r3_t dp_prev; /*previous point*/
    int k;
    bool_t vis_prev = FALSE;   /* Visibility indicator of previous datum in {seq}. */
    bool_t ori_prev = TRUE;    /* Whether previous datum was original or interpolated. */
    start_segment();
    for(k = 0; k < nsmp; k++)
      {
        /* Outputs datum {k} and the rod between datums {k-1} and {k}, if any. */
        /* Get datum {dtk} with index {k} (interpolated): */
        dnae_datum_t *dtk = dnae_seq_get_datum_address(z, k);

        /* Convert to point {dp_this} of {R^3}: */
        r3_t dp_this;
        int j;
        for(j = 0; j < dnae_CHANNELS; j++){
          dp_this.c[j] = dnae_sample_decode(dtk->c[j], z->sfac.f[j]);
          dp_this.c[j] *= magnify;
          dp_this.c[j] += pert->c[j]; 
        }

        /* Attributes of this datum: */
        bool_t vis_this = ((k >= ini) && (k <= fin)); /* Visibility of this datum. */
        bool_t ori_this = (k % nst) == 0; /* Whether this datum is original. */

        if (k > 0)
          { /* End current segment if connecting visible to invisible datums: */
            if ((! vis_this) && vis_prev) { end_segment(); start_segment(); }
            /* Draw rod: */
            extend_segment((vis_this & vis_prev));
            dna_seq_view_tools_draw_stick(&dp_prev, ori_prev, &dp_this, ori_this, cylrad, datrad, quad);
            /* End current segment if connecting from invisible to visible datums: */
            if (vis_this && (! vis_prev)) { end_segment(); start_segment(); }
          }
            
        /* Draw datum: */
        double rad = ( ori_this ? datrad : cylrad );
        /* Original datums are segments by themselves: */
        if (ori_this) { end_segment(); start_segment(); } 
        extend_segment(vis_this); 
        dna_seq_view_tools_draw_ball(&dp_this, rad, quad);
        if (ori_this) { end_segment(); start_segment(); } 
        
        dp_prev = dp_this;
        vis_prev = vis_this;
        ori_prev = ori_this;
      }
    end_segment();
    gluDeleteQuadric(quad);
  }

void dna_seq_view_tools_draw_paired
  ( msm_rung_vec_t *gv,
    dnae_seq_t *x,
    int inix,
    int finx,
    dnae_seq_t *y,
    int iniy,
    int finy,
    double magnify,
    bool_t perturb,
    double radius,
    int nst,
    frgb_t *color_x,
    frgb_t *color_y,
    frgb_t *color_p,
    double dif_cutoff
  )
  {
    double dpert = (perturb ? 0.5*radius : 0.0);
    r3_t x_pert = (r3_t){{ +dpert, +dpert, +dpert }};
    r3_t y_pert = (r3_t){{ -dpert, -dpert, -dpert }};

    /*Draw the sequences*/

    dna_seq_view_tools_draw_sequence(x, magnify, &x_pert, radius, nst, color_x, inix, finx);
    dna_seq_view_tools_draw_sequence(y, magnify, &y_pert, radius, nst, color_y, iniy, finy);

    /*Now we draw the sticks between sequences*/
    dna_seq_view_tools_draw_rungs(x, y, gv, magnify, &x_pert, &y_pert, radius, color_p, dif_cutoff);
  }

void dna_seq_view_tools_draw_rungs
  ( dnae_seq_t *x, 
    dnae_seq_t *y, 
    msm_rung_vec_t *gv, 
    double magnify,
    r3_t *x_pert,
    r3_t *y_pert, 
    double radius,
    frgb_t* color,
    double dif_cutoff
  )
  {
    auto void draw_pairing_stick(int ix, int iy, GLUquadricObj* qd);
    
    void draw_pairing_stick(int ix, int iy, GLUquadricObj* qd)
      { dnae_datum_t dx = dnae_seq_get_datum(x,ix);
        dnae_datum_t dy = dnae_seq_get_datum(y,iy);
        double dif2 = dnae_datum_diffsq(&dx, &(x->sfac), &dy, &(y->sfac));
        assert((dif2 >= 0.0) && (dif2 <= 1.0000000000001));
        double dif = sqrt(dif2);
        r3_t px,py;
        int j;
        for(j = 0; j < 3; j++)
          { px.c[j] = dnae_sample_decode(dx.c[j], x->sfac.f[j]);
            px.c[j] *= magnify;
  	    px.c[j] += x_pert->c[j];
            
	    py.c[j] = dnae_sample_decode(dy.c[j], y->sfac.f[j]);
            py.c[j] *= magnify;
  	    py.c[j] += y_pert->c[j];
          }
        double dim = (dif > dif_cutoff ? 0.5 : 0.0);
        frgb_t gray = (frgb_t){{ 0.500, 0.500, 0.500 }};
        frgb_t col = frgb_mix(dim, &gray, 1-dim, color);
        double rad = (dif > dif_cutoff ? 0.3 : 0.4) * radius;
        glColor3f(col.c[0], col.c[1], col.c[2]);
        dna_seq_view_tools_draw_stick(&px, FALSE, &py, FALSE, rad, 0.0, qd);
      }

    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(color->c[0],color->c[1],color->c[2]);
    GLUquadricObj* quad  = gluNewQuadric();
    msm_rung_t gini = gv->e[0];
    int i;
    for(i = 0; i < gv->ne; i++)
      { msm_rung_t h = gv->e[i];
        int hxk = h.c[0] - gini.c[0];
        int hyk = h.c[1] - gini.c[1];
        draw_pairing_stick(hxk, hyk, quad);
      }
    gluDeleteQuadric(quad);
    glDisable(GL_COLOR_MATERIAL);
  }

void dna_seq_view_tools_draw_ball(r3_t* p, double rad, GLUquadricObj* quad)
  {
    glPushMatrix();
    GLfloat xp = (GLfloat)(p->c[0]);
    GLfloat yp = (GLfloat)(p->c[1]);
    GLfloat zp = (GLfloat)(p->c[2]);
    glTranslatef(xp, yp, zp);
    gluSphere(quad, rad, 6,6);
    glPopMatrix();
  }

void dna_seq_view_tools_draw_stick(r3_t* p, bool_t ptrim, r3_t* q, bool_t qtrim, double rad, double trimrad, GLUquadricObj* quad)
  { 
    /* Get the unit vector from {p} to {q}: */
    r3_t d;
    r3_sub(q, p, &d);
    double dpq = r3_dir(&d,&d);
    
    /* Trim ends is requested: */
    r3_t pp = (*p);
    r3_t qq = (*q);
    if (ptrim) { r3_mix_in(+trimrad, &d, &pp); dpq -= trimrad; }  
    if (qtrim) { r3_mix_in(-trimrad, &d, &qq); dpq -= trimrad; }  
    
    /* Suppress stick if what remained is too short: */
    if (dpq < 1.0e-6) { return; }
  
    glPushMatrix();
    GLfloat xp = (GLfloat)(pp.c[0]);
    GLfloat yp = (GLfloat)(pp.c[1]);
    GLfloat zp = (GLfloat)(pp.c[2]);
    glTranslatef(xp, yp, zp);
    r3_t z = (r3_t){{0,0,1}};
    double angle = (180.0/M_PI) * acos(r3_dot(&z,&d)) ;
    if(angle > 1.0e-3)
      { r3_t t;
        r3_cross(&z, &d, &t);
        GLfloat xt = (GLfloat)(t.c[0]);
        GLfloat yt = (GLfloat)(t.c[1]);
        GLfloat zt = (GLfloat)(t.c[2]);
        glRotatef((GLfloat)angle, xt,yt,zt);
      }
    gluCylinder(quad, rad, rad, (GLfloat)dpq, 6,1);
    glPopMatrix();
  }


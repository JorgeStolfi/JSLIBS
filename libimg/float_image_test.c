/* See {float_image_test.h}. */

/* Last edited on 2024-10-17 15:38:04 by stolfi */ 
/* Created on 2009-06-02 by J. Stolfi, UNICAMP */

#define float_image_test_C_COPYRIGHT \
  "Copyright © 2009  by the State University of Campinas (UNICAMP)"

#include <math.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <r2.h>

#include <float_image.h>
#include <float_image_waves.h>

#include <float_image_test.h>

void float_image_test_gen_ripples(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
  {
    /* Get coordinates of {p}: */
    double x = p->c[0];
    double y = p->c[1];
    
    double srand = M_PI - floor(M_PI); /* Pseudorandom state. */

    auto double frand(void); /* Pseudorandom generator in {[0 _ 1]}. */

    double R = 0.5*hypot(NX, NY); /* Image circumradius. */
    double sigma = 0.30;
    int32_t NW = 5; /* Number of waves. */
    for (uint32_t ic = 0;  ic < NC; ic++) 
      { /* Compute value {fs[ic]} of channel {ic}: */
        double sum_val = 0.0;       /* Sum of values of all component waves. */
        for (uint32_t t = 0;  t < NW; t++)
          { 
            double cx = NX*(0.1 + 0.8*frand()); /* {X} of wave's center. */
            double cy = NY*(0.1 + 0.8*frand()); /* {Y} of wave's center. */
            double rw = R*(0.05 + 0.40*frand()); /* Mid radius of wave packet. */
            double hw = 6.0 + 6.0*frand(); /* Wavelength: */ 
            double amp = 0.7 + 0.3*frand(); /* Max amplitude. */
            
            int32_t T = 1; /* {T=0} open borders, {T=1} torus topology. */
            for (int32_t dx = -T; dx <= +T; dx++)
              for (int32_t dy = -T; dy <= +T; dy++)
                { double r = hypot(x - cx + dx*NX, y - cy + dy*NY); /* Distance from wave center. */
                  double fr = r/rw; /* Distance ralative to packet radius. */
                  if (fr >= 0.01)
                    { double lor = log(fr)/log(2)/sigma;
                      if (fabs(lor) < 8)
                        { double ang = 2*M_PI*(r - rw)/hw; /* Wave argument angle. */
                          double env = exp(-lor*lor/2); /* Packet envelope. */
                          sum_val += amp*env*cos(ang);
                        }
                    }
                }
          }
        /* Map to range {[0 _ 1]}: */
        fs[ic] = (float)(0.5*(erf(sum_val) + 1));
      }
    
    return;
    
    double frand(void)
      { srand = 17*srand;
        srand = srand - floor(srand);
        return srand;
      }
  }
  
void float_image_test_gen_grittie(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
  {
    /* Get coordinates of {p}: */
    double x = p->c[0];
    double y = p->c[1];
  
    int32_t NF = 10;
    double amp[NF];
    double fx[NF];
    double fy[NF]; 
    double phase[NF];
    bool_t verbose = FALSE;
    for (uint32_t ic = 0;  ic < NC; ic++) 
      { 
        float_image_waves_pick(NF, amp, fx, fy, phase, verbose);

        double squash_amp; /* Squash amplitude. */
        /* Choosing a middling squash amplitude: */
        double sum_a2 = 0;
        for (uint32_t kf = 0;  kf < NF; kf++) { sum_a2 += amp[kf]*amp[kf]; }
        double amp_rms = sqrt(sum_a2); /* Root-mean-sum amplitude. */
        squash_amp = 0.5*amp_rms;

        /* Eval the waves: */
        double r = float_image_waves_eval(x, y, NF, amp, fx, fy, phase, squash_amp);

        /* Map {r} to {[0_1]}: */
        r = (r + 1.0)/2;
        fs[ic] = (float)r;
      }
  }  

void float_image_test_gen_stripes(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
  {
    double P = 4.0;  /* Wavelength. */
    for (uint32_t ic = 0;  ic < NC; ic++)
      { int32_t ax = ic % 2; /* Coordinate axis perpendicular to stripes. */
        double wd = (double)(ax == 0 ? NX : NY); /* Image size laong axis {ax}. */
        /* Get coordinate relative to image center: */
        double z = p->c[ax] - wd/2;
        /* Compute the local phase of the wave: */
        double t = 2*M_PI*z/P;
        /* Compute the wave's value: */
        fs[ic] = (float)(0.5*(1 + cos(t)));
      }
  }

void float_image_test_gen_checker(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
  {
    double P0 = 12.0;            /* Wavelength of fundamental. */
    double fr3 = 3.0;            /* Relative frequency of harmonic. */
    double am3 = 0.5/(fr3*fr3);  /* Amplitude of harmonic. */
    
    for (uint32_t ic = 0;  ic < NC; ic++)
      { /* Select the period of the checker depending on the channel: */
        double P = P0*pow(2.0, ic);
        /* Compute the wave's value in the range {[-1 _ +1]}: */
        double v = 1.0;  /* Value of wave. */
        for (uint32_t ax = 0;  ax < 2; ax++)
          { /* Get image size along this axis: */
            double wd = (double)(ax == 0 ? NX : NY);
            /* Get coordinate relative to image center: */
            double z = p->c[ax] - wd/2;
            /* Compute the phase of the fundamental wave (radians): */
            double t = 2*M_PI*z/P;
            /* Compute the wave value: */
            double vax = (1.0 + am3)*cos(t) - am3*cos(fr3*t);
            /* Multiply into {v}: */
            v *= vax;
          }

        /* Adjust {v} to the range {[0 _ 1]}: */
        v = 0.5*(1 + v);
        /* Set into pixel: */
        fs[ic] = (float)v;
      }
  }

void float_image_test_gen_chopsea(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
  {
    /* Get coordinates of {p} relative to image center: */
    double x = p->c[0] - ((double)NX)/2;
    double y = p->c[1] - ((double)NY)/2;
    
    int32_t NW = 45; /* Number of waves. */
    double fCut = 0.25;       /* Cutoff frequency (cycles per pixel). */
    double fBase = 1.173985;  /* Ratio between successive frequencies. */
    for (uint32_t ic = 0;  ic < NC; ic++) 
      { /* Compute value {fs[ic]} of channel {ic}: */
        double sum_val = 0.0;       /* Sum of values of all component waves. */
        double sum_amp = 0.0;       /* Sum of amplitude of all component waves. */
        for (uint32_t t = 0;  t < NW; t++)
          { /* Choose the parameters of wave {t}: */
            double P = 3.0*pow(fBase,t);     /* Wavelength (pixels). */
            int32_t ict = t*NC + ic;         /* A unique integer for channel and component. */
            double azm = 0.618034*ict;       /* Azimuth of wave direction (turns). */ 
            double fx = cos(2*M_PI*azm)/P;   /* X frequency (cycles per pixel). */
            double fy = sin(2*M_PI*azm)/P;   /* Y frequency (cycles per pixel). */
            double phi = 0.618034*(ict+0.5); /* Initial phase of wave (as fraction of cycle). */
            double amp = 1.5 + cos(3.5*ict); /* Unattenuated amplitude. */
            /* Attenuation factor of term {t}: */
            double fr2 = (fx*fx + fy*fy)/(fCut*fCut); /* Freq rel to cutoff, squared. */
            double att = (fr2 >= 1.0 ? 0.0 : 1 - sin(0.5*M_PI*fr2)); /* Attenuation factor. */
            assert(att <= 1.0);
            /* fprintf(stderr, "    %3d  fr = %10.5f %10.5f", t,  fx, fy); */
            /* fprintf(stderr, "  P = %10.5f", P); */
            /* fprintf(stderr, "  phi = %10.5f  amp = %10.5f att = %10.5f\n", phi, amp, att); */

            /* Wave phase of term {t} at {x,y} (radians). */
            double tph = 2*M_PI*(fx*x + fy*y + phi); 
            /* Value of term {t}: */
            double val = att*amp*sin(tph); 
            /* Accumulate: */
            sum_val += amp*val;
            sum_amp += amp;
          }
        /* Map to range {[0 _ 1]}: */
        fs[ic] = (float)((sum_val/sum_amp + 1)/2);
      }
  }

void float_image_test_gen_bullsex(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
  {
    /* Get coordinates of {p} relative to image center: */
    double x = p->c[0] - ((double)NX)/2;
    double y = p->c[1] - ((double)NY)/2;
    
    /* Adjust parameters {A} and {B} so that the wavelength is {L} pixels at farthest edge: */
    double L = 4.0;
    double H = fmax(NX, NY)/2; 
    double R = hypot(NX, NY)/2;
    double s = R/H;
    double B = log(s)/(R - H);
    double A = M_PI/(B*sinh(B*H))/L;
    
    double r = hypot(x, y);
    double arg = A*(cosh(B*r) - 1);
    for (uint32_t ic = 0;  ic < NC; ic++)
      { /* Vary the phase of the wave according to the channel: */
        double phase = (0.25 + 0.50*ic)*M_PI;
        fs[ic] = (float)(0.5*(cos(arg + phase) + 1)); 
      }
  }

void float_image_test_gen_bullsqr(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
  {
    /* Get coordinates of {p} relative to image center: */
    double x = p->c[0] - ((double)NX)/2;
    double y = p->c[1] - ((double)NY)/2;
    
    /* Determine {A} so that the wavelength is {L} pixels at the edge: */
    double L = 4.0;
    double H = fmax(NX, NY)/2;
    double A = M_PI/(2*H)/L;
    
    double r = hypot(x, y);
    double arg = A*r*r;
    for (uint32_t ic = 0;  ic < NC; ic++)
      { /* Vary the phase of the wave according to the channel: */
        double phase = (0.25 + 0.50*ic)*M_PI;
        fs[ic] = (float)(0.5*(cos(arg + phase) + 1)); 
      }
  }

void float_image_test_paint
  ( float_image_t *img, 
    float_image_test_generator_t *proc,
    int32_t ns
  )
  {
    /* Get image dimensions: */
    int32_t NC = (int32_t)img->sz[0];
    int32_t NX = (int32_t)img->sz[1];
    int32_t NY = (int32_t)img->sz[2];
    
    float v[NC];  /* Image pixel samples: */
    double vd[NC];  /* Accumulated subsamples in a pixel: */

    /* Scan pixels: */
    double ns2 = ns*ns; /* Total subsamples in each pixel. */
    for (uint32_t iy = 0;  iy < NY; iy++)
      { for (uint32_t ix = 0;  ix < NX; ix++)
          { /* Accumulated pixel value: */
            for (uint32_t ic = 0;  ic < NC; ic++) { vd[ic] = 0.0; }
            /* Evaluate procedure at subsampling points and accumulate: */
            for (uint32_t iys = 0;  iys < ns; iys++)
              { for (uint32_t ixs = 0;  ixs < ns; ixs++)
                  { /* Get coordinates of sampling point in pixel, relative to image center: */
                    r2_t p = (r2_t){{ ix + (ixs + 0.5)/ns, iy + (iys + 0.5)/ns }};
                    /* Evaluate procedural image: */
                    proc(&p, NC, NX, NY, v);
                    /* Accumulate sample values: */
                    for (uint32_t ic = 0;  ic < NC; ic++) { vd[ic] += (double)v[ic]; }
                  }
              }
            /* Reduce sum to average and store in pixel: */
            for (uint32_t ic = 0;  ic < NC; ic++) { v[ic] = (float)(vd[ic]/ns2); }
            float_image_set_pixel(img, ix, iy, v);
          }
      }
  }

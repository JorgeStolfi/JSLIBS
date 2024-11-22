/* Last edited on 2024-11-20 13:23:45 by stolfi */
/* Extra tests when these functions are implemented. */

void test_hr3_pmap_congruence_from_point_and_dir(bool_t verbose, bool_t flip);
void test_hr3_pmap_similarity_from_two_points(bool_t verbose, bool_t flip);
void test_hr3_pmap_rotation_and_scaling(bool_t verbose);

void test_hr3_pmap_diff_sqr(bool_t verbose);
void test_hr3_pmap_mismatch_sqr(bool_t verbose);
void test_hr3_pmap_aff_discr_sqr(bool_t verbose);
void test_hr3_pmap_deform_sqr(bool_t verbose);

    test_hr3_pmap_diff_sqr(verbose);
    test_hr3_pmap_deform_sqr(verbose);
    test_hr3_pmap_mismatch_sqr(verbose);

    test_hr3_pmap_aff_discr_sqr(verbose);

    test_hr3_pmap_congruence_from_point_and_dir(verbose, FALSE);
    test_hr3_pmap_congruence_from_point_and_dir(verbose, TRUE);
    test_hr3_pmap_similarity_from_two_points(verbose, FALSE);
    test_hr3_pmap_similarity_from_two_points(verbose, TRUE);
    test_hr3_pmap_rotation_and_scaling(verbose);
 
void test_hr3_pmap_diff_sqr(bool_t verbose)
  { if (verbose) { fprintf(stderr, "--- hr3_pmap_diff_sqr ---\n"); }
    hr3_pmap_t M, N;
    /* Throw two maps and normalizes their matrices: */
    for (int32_t k = 0; k < 2; k++)
      { N = M;
        h2tt_throw_pmap(&M); 
        r4x4_normalize(&(M.dir));
        r4x4_normalize(&(M.inv));
      }
    /* Compute exepected diff sqr {d2_exp}: */
    double d2_exp = 0;
    for (int32_t i = 0; i < NH; i++)
      { for (int32_t j = 0; j < NH; j++)
          { double d = M.dir.c[i][j] - N.dir.c[i][j];
            d2_exp += d*d;
            d = M.inv.c[i][j] - N.inv.c[i][j];
            d2_exp += d*d;
          }
      }
    /* Scale the matrices by arbitrary amounts: */
    for (int32_t i = 0; i < NH; i++)
      { for (int32_t j = 0; j < NH; j++)
          { M.dir.c[i][j] *= 2;
            M.inv.c[i][j] *= 4;
            N.inv.c[i][j] *= 0.25;
            N.dir.c[i][j] *= 0.0625;
          }
      }
    /* Compute {hr3_pmap_diff_sqr} and compare: */
    double d2_cmp = hr3_pmap_diff_sqr(&M, &N);
    h2tt_check_eps(d2_cmp, d2_exp, 1.0e-13, "hr3_pmap_diff_sqr failed");
  }

void test_hr3_pmap_deform_sqr(bool_t verbose)
  { 
    fprintf(stderr, "!! {hr3_pmap_deform_sqr} NOT TESTED\n");
  }

void test_hr3_pmap_mismatch_sqr(bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_pmap_mismatch_sqr ---\n"); }
    
    bool_t debug = FALSE;

    /* Choose a bunch of points: */
    uint32_t np = 7;
    r3_t p1[np], p2[np];
    for (int32_t ip = 0; ip < np; ip++)
      { r3_throw_cube(&(p1[np]));
        r3_throw_cube(&(p2[np]));
      }
    /* Choose a projective map: */
    hr3_pmap_t M; h2tt_throw_pmap(&M);
    
    if (debug) { hr3_pmap_print(stderr, &M, "  M =\n", "\n"); }
    
    /* Compute the expected mismatch {m2_exp}: */
    double sum2 = 0.0;
    hr3_pmap_t N = hr3_pmap_inv(&M);
    for (int32_t ip = 0; ip < np; ip++)
      { r3_t *p1k = &(p1[ip]);
        r3_t *p2k = &(p2[ip]);
        r3_t q1k = hr3_pmap_r3_point(p1k, &M);
        r3_t q2k = hr3_pmap_r3_point(p2k, &N);
        double d2 = r3_dist_sqr(&q1k, &q2k);
        if (debug)
          { r3_gen_print(stderr, p1k, "%+8.5f", "  ( ", " ", " )");
            r3_gen_print(stderr, &q1k, "%+12.8f", " -> ( ", " ", " )");
            fprintf(stderr, " |%.6f| ", d2);
            r3_gen_print(stderr, &q2k, "%+12.8f", "( ", " ", " ) <- ");
            r3_gen_print(stderr, p2k, "%+8.5f", "( ", " ", " )\n");
          }
        sum2 += d2;
      }
    double m2_exp = sum2/np;
    
    /* Compute the mismatch by the library: */
    double m2_cmp = hr3_pmap_mismatch_sqr(&M, np, p1, p2);
    h2tt_check_eps(m2_cmp, m2_exp, 1.0e-13, "hr3_pmap_diff_sqr");
  }

void test_hr3_pmap_aff_discr_sqr(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_aff_discr_sqr ---\n"); }
    
    bool_t debug = FALSE;

    hr3_pmap_t M, N;
     h2tt_throw_aff_map(&M);
     h2tt_throw_aff_map(&N);
    if (debug)
      { h2tt_print_pmap("M", &M);
        h2tt_print_pmap("N", &N);
      }
    double mis2 = hr3_pmap_aff_discr_sqr(&M, &N);
    /* The integrand |(M - N)(u(t))|^2 is a sinusoidal function with freqs 0,1,2. */
    /* Can be integrated numerically with 5 or more samples. */
    uint32_t nang = 7;
    double sum_d2 = 0; 
    for (int32_t i = 0; i < nang; i++)
      { double ang = 2*M_PI*((double)i)/((double)nang);
        r3_t p = (r3_t){{ cos(ang), sin(ang) }};
        r3_t q = hr3_pmap_r3_point(&p, &M);
        r3_t r = hr3_pmap_r3_point(&p, &N);
        double d2 = r3_dist_sqr(&q, &r);
        sum_d2 += d2;
      }
    double cis2 = sum_d2/nang;
    h2tt_check_num_eps("mis", mis2, cis2, 0.0000001, "hr3_pmap_aff_discr_sqr failed");
  }

void test_hr3_pmap_congruence_from_point_and_dir(bool_t verbose, bool_t flip)
  { 
    if (verbose) { fprintf(stderr, "--- hr3_pmap_congruence_from_point_and_dir ---\n"); }

    bool_t debug = FALSE;
        
    r3_t oc = (r3_t){{ 0.0, 0.0 }};
    r3_t pc = (r3_t){{ 1.0, 0.0 }};
    r3_t qc = (r3_t){{ 0.0, 1.0 }};
    r3_t rc = (r3_t){{ 1.0, 1.0 }};

    r3_t ocM; r3_throw_cube(&ocM);
    r3_t udM; r3_throw_dir(&udM);

    /* Determine the vector {v} that is going to be the image of vector {(0,1)}: */
    r3_t vdM = (r3_t){{ -udM.c[1], +udM.c[0] }}; 
    if (flip) { r3_neg(&vdM, &vdM); }
    
    if (debug) 
      { fprintf(stderr, "  flip = %c\n", "FT"[flip]);
        r3_gen_print(stderr, &ocM, "%12.8f", "  ocM = ( ", " ", " )\n");
        r3_gen_print(stderr, &udM, "%12.8f", "  udM = ( ", " ", " )\n");
        r3_gen_print(stderr, &vdM, "%12.8f", "  vdM = ( ", " ", " )\n");
      }
    
    /* Compute the expected images of points {pc,qc,rc}: */
    r3_t pcM; r3_add(&ocM, &udM, &pcM);
    r3_t qcM; r3_add(&ocM, &vdM, &qcM);
    r3_t rcM; r3_add(&pcM, &vdM, &rcM);

    hr3_pmap_t A = hr3_pmap_congruence_from_point_and_dir(&ocM, &udM, flip);

    /* Check whether the map works: */
    h2tt_check_pmap_r3_point("o", &oc, &A, FALSE, &ocM, "hr3_pmap_congruence_from_point_and_dir failed");
    h2tt_check_pmap_r3_point("p", &pc, &A, FALSE, &pcM, "hr3_pmap_congruence_from_point_and_dir failed");
    h2tt_check_pmap_r3_point("q", &qc, &A, FALSE, &qcM, "hr3_pmap_congruence_from_point_and_dir failed");
    h2tt_check_pmap_r3_point("r", &rc, &A, FALSE, &rcM, "hr3_pmap_congruence_from_point_and_dir failed");
  }
    
void test_hr3_pmap_similarity_from_two_points(bool_t verbose, bool_t flip)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_similarity_from_two_points ---\n"); }

    bool_t debug = FALSE;
        
    r3_t oc = (r3_t){{ 0.0, 0.0 }};
    r3_t pc = (r3_t){{ 1.0, 0.0 }};
    r3_t qc = (r3_t){{ 0.0, 1.0 }};
    r3_t rc = (r3_t){{ 1.0, 1.0 }};

    r3_t ocM; r3_throw_cube(&ocM);
    r3_t pcM; r3_throw_cube(&pcM);
        
    /* Determine the vector {vdM} that is going to be the image of vector {(0,1)}: */
    r3_t udM; r3_sub(&pcM, &ocM, &udM);
    r3_t vdM = (r3_t){{ -udM.c[1], +udM.c[0] }};
    if (flip) { r3_neg(&vdM, &vdM); }
    
    if (debug) 
      { fprintf(stderr, "  flip = %c\n", "FT"[flip]);
        r3_gen_print(stderr, &ocM, "%12.8f", " ocM = ( ", " ", " )\n");
        r3_gen_print(stderr, &pcM, "%12.8f", " pcM = ( ", " ", " )\n");
        r3_gen_print(stderr, &udM, "%12.8f", " udM =  ( ", " ", " )\n");
        r3_gen_print(stderr, &vdM, "%12.8f", " vdM =  ( ", " ", " )\n");
      }
    
    /* Compute the expected images of points {pc,qc,rc}: */
    r3_t qcM; r3_add(&ocM, &vdM, &qcM);
    r3_t rcM; r3_add(&pcM, &vdM, &rcM);
    
    hr3_pmap_t A = hr3_pmap_similarity_from_two_points(&ocM, &pcM, flip);

    /* Check whether the map works: */
    h2tt_check_pmap_r3_point("o", &oc, &A, FALSE, &ocM, "hr3_pmap_similarity_from_two_points failed");
    h2tt_check_pmap_r3_point("p", &pc, &A, FALSE, &pcM, "hr3_pmap_similarity_from_two_points failed");
    h2tt_check_pmap_r3_point("q", &qc, &A, FALSE, &qcM, "hr3_pmap_similarity_from_two_points failed");
    h2tt_check_pmap_r3_point("r", &rc, &A, FALSE, &rcM, "hr3_pmap_similarity_from_two_points failed");
  }


void test_hr3_pmap_rotation_and_scaling(bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- hr3_pmap_rotation_and_scaling ---\n"); }
    double ang = dabrandom(-1.58, +1.58);
    double scale = dabrandom(0.5, 2.0);
    hr3_pmap_t M = hr3_pmap_rotation_and_scaling(ang, scale);
    double ca = cos(ang);
    double sa = sin(ang);
    for (int32_t k = 0; k < 5; k++)
      { r3_t p; r3_throw_cube(&p);
        r3_t q; 
        q.c[0] = + ca*scale*p.c[0] - sa*scale*p.c[1];
        q.c[1] = + sa*scale*p.c[0] + ca*scale*p.c[1];
        h3tt_check_pmap_r3_point("p", &p, &M, FALSE, &q, "hr3_pmap_rotation_and_scaling failed");
      }
  }

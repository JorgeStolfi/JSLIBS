#include <stdint.h>

void foo (
    Float a_ctr
,    Float a_0
,    Float b_ctr
,    Float b_1
,    Float c_ctr
,    Float c_2
,    Float *ra_ctr,
    Float *ra_3
,    int32_t *nops
  )
  {
    Float alpha, beta, zeta, gamma, delta, tmp;
    Float d_ctr, s_rho;    Interval d_range;    Float d_3;
    Float e_ctr, s_rho;    Interval e_range;    Float f_ctr, s_rho;    Interval f_range;    Float f_4;
    Float g_ctr, s_rho;    Interval g_range;    Float g_4;
    Float g_5;
    Float h_ctr, s_rho;    Interval h_range;    Float h_3;
    Float i_ctr, s_rho;    Interval i_range;    Float i_3;
    Float i_4;
    Float i_5;
    Float i_6;
    d_3 = Zero;
    d_rho = Zero;
    flt_add(a_ctr, b_ctr, &d_ctr, &d_3);
    flt_add(a_0, b_0, &tmp, &d_3);
    ROUND_UP; d_3 += FABS(tmp);
    flt_add(a_1, b_1, &tmp, &d_3);
    ROUND_UP; d_3 += FABS(tmp);
    d_rho += FABS(d_3);
    d_range = ia_meet(ia_add(a_range, b_range), ia_const(d_ctr, d_rho))
    e_ctr = (Float) 3;
    e_rho = Zero;
    e_range = (Interval) {e_ctr, e_ctr};
    cheb_sqrt(c_range, &alpha, &zeta, &gamma, &delta);
    f_4 = delta;
    flt_affine(c_ctr, alpha, zeta, gamma, &f_ctr, &f_4);
    flt_scale(c_2, alpha, zeta, &tmp, &f_4);
    f_4 += FABS(tmp);
    f_range = ia_sqrt(c_range);
    g_5 = Zero;
    g_rho = Zero;
    flt_add(e_ctr, f_ctr, &g_ctr, &g_5);
    flt_add(e_4, f_4, &g_4, &g_5);
    ROUND_UP; g_rho += FABS(g_4);
    g_rho += FABS(g_5);
    g_range = ia_meet(ia_add(e_range, f_range), ia_const(g_ctr, g_rho))
    h_ctr = - d_ctr;
    h_3 = - d_3;
    h_range = ia_neg(d_range);
    i_6 = Zero;
    i_rho = Zero;
    flt_add(g_ctr, h_ctr, &i_ctr, &i_6);
    flt_add(g_3, h_3, &i_3, &i_6);
    ROUND_UP; i_rho += FABS(i_3);
    flt_add(g_4, h_4, &i_4, &i_6);
    ROUND_UP; i_rho += FABS(i_4);
    flt_add(g_5, h_5, &i_5, &i_6);
    ROUND_UP; i_rho += FABS(i_5);
    i_rho += FABS(i_6);
    i_range = ia_meet(ia_add(g_range, h_range), ia_const(i_ctr, i_rho))
    ROUND_UP;
    *ra_ctr = i_ctr;
    *ra_9 = Zero;
    *ra_9 += ABS(i_3);
    *ra_9 += ABS(i_4);
    *ra_9 += ABS(i_5);
    *ra_9 += ABS(i_6);
  }


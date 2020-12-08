/* See {nmsim_class_neuron.h} */
/* Last edited on 2020-12-04 21:01:05 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <vec.h>
#include <filefmt.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <nmsim_basic.h>

#include <nmsim_firing_func.h>
#include <nmsim_class_neuron.h>
 
nmsim_class_neuron_t* nmsim_class_neuron_new
  ( double V_B,   
    double V_R,   
    double c_B,   
    double M_R,   
    double M_mu,  
    double H_R,   
    double H_mu,  
    struct nmsim_firing_func_t *Phi
  )
  {
    nmsim_class_neuron_t *nclass = notnull(malloc(sizeof(nmsim_class_neuron_t)), "no mem");
    (*nclass) = (nmsim_class_neuron_t)
      { .V_B = V_B, .V_R = V_R, 
        .c_B = c_B, 
        .M_R = M_R, .M_mu = M_mu,
        .H_R = H_R, .H_mu = H_mu,
        .Phi = (*Phi)
      };
    return nclass;
  }

void nmsim_class_neuron_free(nmsim_class_neuron_t *nclass)
  { 
    if (nclass != NULL) { free(nclass); }
  }

double nmsim_class_neuron_compute_M(nmsim_class_neuron_t *nclass, nmsim_step_count_t age)
  {
    demand(age >= 0, "invalid firing age");
    double M = nclass->M_R;
    if ((age > 0) && (M != 1.0))
      { double mu = nclass->M_mu;
        if (mu == 0.0)
          { M = 1.0; }
        else if (mu != 1.0)
          { M = 1.0 - (1.0 - M)*pow(mu, (double)age); }
      }
    return M;
  }

double nmsim_class_neuron_compute_H(nmsim_class_neuron_t *nclass, nmsim_step_count_t age)
  {
    demand(age >= 0, "invalid firing age");
    double H = nclass->H_R;
    if ((age > 0) && (H != 1.0))
      { double mu = nclass->H_mu;
        if (mu == 0.0)
          { H = 1.0; }
        else if (mu != 1.0)
          { H = 1.0 - (1.0 - H)*pow(mu, (double)age); }
      }
    return H;
  }
    
double nmsim_class_neuron_recharge(nmsim_class_neuron_t *nclass, double V, double M)
  {
    double V_B = nclass->V_B;
    double c = nclass->c_B * M;
    return V_B + (V - V_B)*c;
  }
  
void nmsim_class_neuron_throw_state
  ( nmsim_class_neuron_t *nclass, 
    double *VP, 
    nmsim_step_count_t *ageP
  )
  { /* How could it be improved? */
    (*VP) = 0.5*(nclass->V_R + nclass->V_B);
    (*ageP) = (nmsim_step_count_t)int32_abrandom(0, 9);
  }

void nmsim_class_neuron_write(FILE *wr, nmsim_class_neuron_t *nclass, double timeStep)
  {
    char *ind1 = "      ";   /* Indent of header lines. */
    char *ind2 = "        "; /* Indent of parameter lines. */
    
    auto void write_tau_param_from_mu(char *name, double mu);
      /* Writes a line "{name} = {tau}", where {tau} is the
        characteristic decay time that corresponds to the
        decay fator {mu} in a simulation with the given {timeStep}. */
    
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_class_neuron_FILE_TYPE, nmsim_class_neuron_VERSION);
    nmsim_write_double_param      (wr, ind2, "V_R",  nclass->V_R, nmsim_write_VIJ_PREC, TRUE,  TRUE,  FALSE);
    nmsim_write_double_param      (wr, ind2, "V_B",  nclass->V_B, nmsim_write_VIJ_PREC, TRUE,  TRUE,  FALSE);
    write_tau_param_from_mu ("V_tau",  nclass->c_B);
    nmsim_write_double_param      (wr, ind2, "M_R",  nclass->M_R, nmsim_write_MH_PREC,  FALSE, TRUE,  TRUE );
    write_tau_param_from_mu ("M_tau", nclass->M_mu);
    nmsim_write_double_param      (wr, ind2, "H_R",  nclass->H_R, nmsim_write_MH_PREC,  FALSE, TRUE,  TRUE );
    write_tau_param_from_mu ("H_tau", nclass->H_mu);
    
    /* Write the firing function parameters: */
    fprintf(wr, "%sPhi_class = %c\n", ind2, nclass->Phi.class);
    nmsim_write_double_param(wr, ind2, "V_M",  nclass->Phi.V_M, nmsim_write_VIJ_PREC, TRUE,  TRUE,  FALSE);
    nmsim_write_double_param(wr, ind2, "V_D",  nclass->Phi.V_D, nmsim_write_VIJ_PREC, FALSE, TRUE,  FALSE);
    
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_class_neuron_FILE_TYPE);

    fflush(wr);
    
    return;
    
    /* INTERNAL IMPLEMENTATIONS */

    void write_tau_param_from_mu(char *name, double mu)
      { double tau;
        if (mu == 0)
          { tau = 0; }
        else if (mu == 1.0)
          { tau = INF; }
        else
         { tau = -timeStep/log(mu); }
        nmsim_write_double_param(wr, ind2, name, tau, nmsim_write_tau_PREC, FALSE, TRUE, FALSE);
      }
  }

nmsim_class_neuron_t *nmsim_class_neuron_read(FILE *rd, double timeStep)
  {
    auto double read_mu_from_tau_param(char *name);
      /* Reads a line \"{name} = {value}\", including the end-of-line,
       and then converts {value} from a characteristic time
       to a decay factor based on the given {timeStep}. */
    
    /* Read header line: */
    filefmt_read_header(rd, nmsim_class_neuron_FILE_TYPE, nmsim_class_neuron_VERSION);
    
    /* Read the parameters: */
    double V_R  = nmsim_read_double_param(rd, "V_R", -200.0, +200.0);
    double V_B  = nmsim_read_double_param(rd, "V_B", -200.0, +200.0);
    double c_B  = read_mu_from_tau_param("V_tau");
    double M_R  = nmsim_read_double_param(rd, "M_R", 0.0, 10.0);
    double M_mu = read_mu_from_tau_param("M_tau");
    double H_R  = nmsim_read_double_param(rd, "H_R", 0.0, 10.0);
    double H_mu = read_mu_from_tau_param("H_tau");
    
    /* Read name and parameters of the firing function: */
    nmsim_firing_func_class_t Phi_class = nget_char(rd, "Phi_class"); fget_eol(rd);
    double V_M = nmsim_read_double_param(rd, "V_M", -200.0, +200.0);
    double V_D = nmsim_read_double_param(rd, "V_D", 0.0, +200.0);
    nmsim_firing_func_t Phi = nmsim_firing_func_make(Phi_class, V_M, V_D);
    
    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_class_neuron_FILE_TYPE);
    
    /* Build record: */
    nmsim_class_neuron_t *nclass = nmsim_class_neuron_new
      ( V_B, V_R, c_B, M_R, M_mu, H_R, H_mu, &Phi );

    return nclass;

    /* INTERNAL IMPLEMENTATIONS */

    double read_mu_from_tau_param(char *name)
      { double tau = nmsim_read_double_param(rd, name, 0.0, INF);
        if (tau < 0.02*timeStep)
          { return 0.0; }
        else if (tau == INF)
          { return 1.0; }
        else
          { return exp(-timeStep/tau); } 
      }
  }

void nmsim_class_neuron_show(FILE *wr, char *pref, nmsim_class_neuron_t *nclass, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf(wr, "{ VB = %+.2f VR = %+.2f", nclass->V_B, nclass->V_R);
    fprintf(wr, " cB = %.6f", nclass->c_B);
    if (nclass->M_R != 1.0) 
      { fprintf(wr, " MR = %.6f Mmu = %.6f", nclass->M_R, nclass->M_mu); }
    if (nclass->H_R != 1.0)
      { fprintf(wr, " HR = %.6f Hmu = %.6f", nclass->H_R, nclass->H_mu); }
    fprintf(wr, " Phi = %c", nclass->Phi.class);
    fprintf(wr, " VM = %+.2f VD = %.2f", nclass->Phi.V_M, nclass->Phi.V_D);
    fprintf(wr, " }");
    if (suff != NULL) { fputs(suff, wr); }
  }

nmsim_class_neuron_t* nmsim_class_neuron_throw(void)
  { 
    double V_R = nmsim_throw_double(-70.0, -50.0);
    double V_B = nmsim_throw_double(-80.0, -60.0);
    double c_B = nmsim_throw_double(0.750, 0.950);
    double M_R = nmsim_throw_double(0.300, 0.800);
    double M_mu = nmsim_throw_double(0.750, 0.950);
    double H_R = nmsim_throw_double(0.300, 0.800);
    double H_mu = nmsim_throw_double(0.750, 0.950);
    nmsim_firing_func_class_t Phi_class = "GLN"[int32_abrandom(0,1)];
    double V_M = nmsim_throw_double(-50.0, -30.0);
    double V_D = nmsim_throw_double(5.0, 25.0);
    
    nmsim_firing_func_t Phi = nmsim_firing_func_make(Phi_class, V_M, V_D);
    
    nmsim_class_neuron_t *nclass = nmsim_class_neuron_new
      ( V_R, V_B, c_B, M_R, M_mu, H_R, H_mu, &Phi);
      
    return nclass;
  }

void nmsim_class_neuron_compare(nmsim_class_neuron_t *nclass_read, nmsim_class_neuron_t *nclass_orig)
  {
    nmsim_compare_double_param("V_R",  nclass_read->V_R,  nclass_orig->V_R,  nmsim_write_VIJ_PREC, TRUE, FALSE);
    nmsim_compare_double_param("V_B",  nclass_read->V_B,  nclass_orig->V_B,  nmsim_write_VIJ_PREC, TRUE, FALSE);
    nmsim_compare_double_param("c_B",  nclass_read->c_B,  nclass_orig->c_B,  nmsim_write_MH_PREC,  TRUE, TRUE );
    nmsim_compare_double_param("M_R",  nclass_read->M_R,  nclass_orig->M_R,  nmsim_write_MH_PREC,  TRUE, TRUE );
    nmsim_compare_double_param("M_mu", nclass_read->M_mu, nclass_orig->M_mu, nmsim_write_MH_PREC,  TRUE, TRUE );
    nmsim_compare_double_param("H_R",  nclass_read->H_R,  nclass_orig->H_R,  nmsim_write_MH_PREC,  TRUE, TRUE );
    nmsim_compare_double_param("H_mu", nclass_read->H_mu, nclass_orig->H_mu, nmsim_write_MH_PREC,  TRUE, TRUE );

    nmsim_firing_func_compare(&(nclass_read->Phi), &(nclass_orig->Phi));
  }    

vec_typeimpl(nmsim_class_neuron_ref_vec_t,nmsim_class_neuron_ref_vec,nmsim_class_neuron_ref_t);

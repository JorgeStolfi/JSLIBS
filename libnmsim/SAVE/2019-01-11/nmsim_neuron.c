/* See {nmsim_neuron_elem.h} */
/* Last edited on 2016-08-11 11:00:21 by stolfilocal */

#define _GNU_SOURCE

#include <nmsim_firing_func.h>
#include <nmsim_neuron_elem.h>
 
nmsim_neuron_class_t* nmsim_neuron_class_new
  ( double V_R,   /* Potential after firing. */
    double mu_V,  /* Potential decay factor. */
    double G_R,   /* Firing gain factor after firing. */
    double mu_G,  /* Firing gain recovery factor. */
    double H_R,   /* Output gain factor after firing. */
    double mu_H,  /* Output gain recovery factor. */
    struct nmsim_firing_func_t *Phi
  )
  {
    nmsim_neuron_class_t* pm = notnull(malloc(sizeof(nmsim_neuron_class_t)), "no mem");
    (*pm) = (nmsim_neuron_class_t)
      { .V_R = V_R, .mu_V = mu_V,
        .G_R = G_R, .mu_G = mu_G,
        .H_R = H_R, .mu_H = mu_H,
        .Phi = Phi
      };
    return pm;
  }

nmsim_neuron_class_read(FILE *rd, double timeStep)
  {
    auto double read_param(char *name, double vmin, double vmax);
      /* Reads a line \"{name} = {value}\", including the end-of-line,
        and checks that {value} is in the range {[vmin_vmax]}. */

    auto double read_mu_from_tau_param(char *name, double vmin, double vmax);
      /* Reads a line \"{name} = {value}\", including the end-of-line,
       checks that {value} is in the range {[vmin_vmax]},
       and then converts {value} from a chaacteristic time
       to a decay factor based on the given {timeStep}. */
    
    double read_param(char *name, double vmin, double vmax)
      { double v = nget_double(rd, name);
        if ((v < vmin) || (v > vmax))
          { fprintf(stderr, "** parameter {%s} = %24.16e is out of range [ %24.16e _ %24.16e ]\n", name, v, vmin, vmax); 
            demand(FALSE, "aborted");
          }
        fget_eol(rd);
        return v;
      }
    
    double read_mu_from_tau_param(char *name)
      { double tau = read_param(name, 0.0, INF);
        if (tau < 0.02*timeStep)
          { return 0.0; }
        else if (tau == INF)
          { return 1.0; }
        else
          { return exp(-timeStep/tau); } 
      }

    /* Read header line: */
    filefmt_read_header(rd, "nmsim_neuron_elem_parms", nmsim_neuron_class_VERSION);
    
    /* Read the parameters: */
    double V_R  = read_param("V_R", -200.0, +200.0, );
    double mu_V = read_mu_from_tau_param("tau_V");
    double G_R  = read_param("G_R", 0.0, 1.0);
    double mu_G = read_mu_from_tau_param("tau_G");
    double H_R  = read_param("H_R", 0.0, 1.0);
    double mu_H = read_mu_from_tau_param("tau_H");
    
    /* Read name and parameters of the firing function: */
    
    char *Phi_name
    


    /* Read footer line: */
    filefmt_read_footer(rd, type);
  }

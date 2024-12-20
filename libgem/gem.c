/* See {gem.h} */
/* Last edited on 2024-12-05 10:25:50 by stolfi */

#define gem_C_copyright "Copyright © 2014 State University of Campinas (UNICAMP)"

#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <gem.h>

typedef struct gem_node_t
  { int dim;           /* Allocated dimension, in {0..gem_MAX_DIM}. */
    int data;          /* Used to associate data to the node; {-1} if no data. */
    int label;         /* Used in gem traversal. */
    gem_ref_t adj[1];  /* Must be the last field. Actually {.dim+1} elements. */
  } gem_node_t;

bool_t gem_node_is_marked(gem_ref_t p, int nn, gem_ref_t node[]);
  /* Used by {gem_domain_traverse} to determine if a node {p} has been visited
    already or is in the queue {node[0..nn-1]}, waiting to be visited.  Specifically,
    if the label {k} of {p} is in the range {0..nn-1}, and {p == node[k]}. */

gem_ref_t gem_node_new(int d)
  { demand((d >= 0) && (d <= gem_DIM_MAX), "invalid dimension");
    /* Allocate {d+1} adjacency links (one in the record type, plus {d}): */
    gem_ref_t p = notnull(malloc(sizeof(gem_node_t) + d*sizeof(gem_ref_t)), "no mem");
    p->dim = d;
    p->data = -1;
    p->label = 0;
    /* Set all links to itself (all walls are on the free border): */
    int i;
    for (i = 0; i <= d; i++) { p->adj[i] = p; }
    return p;
  }

int gem_node_dim(gem_ref_t p)
  { return (p == NULL ? -1 : p->dim); }

void gem_node_free(gem_ref_t p)
  { if (p == NULL) { return; }
    /* Detach the node: */
    int d = p->dim;
    int i;
    for (i = 0; i <= d; i++) { gem_splice(p, p->adj[i], i); }
    free(p);
  }

#define gem_INI_NODES 500  
  /* Initial size of node stacks. */
    
void gem_component_free(gem_ref_t p)
  { if (p == NULL) { return; }
    gem_ref_vec_t vis = gem_ref_vec_new(gem_INI_NODES);
    int nvis = 0;
    gem_traverse(p, gem_DIM_MAX, &vis, &nvis);
    int i;
    for (i = 0; i < nvis; i++) { gem_node_free(vis.e[i]); }
    free(vis.e);
  }

gem_ref_t gem_step(gem_ref_t p, int i)
  { demand(i >= 0, "invalid link index");
    return (i > p->dim ? p : p->adj[i]);
  }

void gem_set_data(gem_ref_t p, int data)
  { p->data = data; }

int gem_get_data(gem_ref_t p)
  { return p->data; }

int gem_get_label(gem_ref_t p)
  { return p->label; }

void gem_splice(gem_ref_t a, gem_ref_t b, int i)
  { demand((a != NULL) && (b != NULL), "null node pointer");
    demand((i >= 0) && (i <= a->dim) && (i <= b->dim), "invalid link index");
    gem_ref_t a1 = a->adj[i];
    gem_ref_t b1 = b->adj[i];
    demand(((a1 == a) && (b1 == b)) || ((a1 == b) && (b1 == a)), "invalid call");
    a->adj[i] = b1;
    b->adj[i] = a1;
  }

bool_t gem_node_is_marked(gem_ref_t p, int nn, gem_ref_t node[])
  {
    return ((p->label >= 0) && (p->label < nn) && (node[p->label] == p));
  }

void gem_domain_traverse(gem_ref_t root, int nr, int R[], gem_ref_vec_t *visP, int *nvisP)
  { demand(root != NULL, "null {root}");
    demand(visP != NULL, "null {visP}");
    demand(nvisP != NULL, "null {nvisP}"); 
    int nvis = (*nvisP);  /* Number of nodes seen so far (completely or partially). */
    demand((nvis >= 0) && (nvis <= visP->ne), "invalid {*nvisP}");
    
    if (! gem_node_is_marked(root, nvis, visP->e))
      { int tvis = nvis;/* Number of nodes whose R-links have been followed. */
        /* Stack the root: */
        gem_ref_vec_expand(visP, nvis);
        visP->e[nvis] = root; root->label = nvis;  nvis++;
        /* Propagate: */
        while (tvis < nvis)
          { /* Grab the next incompletely visited nodeP {p}: */
            gem_ref_t p = visP->e[tvis];;
            /* Check its R-links: */
            int jr;
            for (jr = 0; jr < nr; jr++) 
              { /* Get the {R[jr]}-neighbor {q} of {p}, stack it if not marked: */
                int i = (R == NULL ? jr : R[jr]);
                if (i <= p->dim)
                  { gem_ref_t q =  p->adj[i];
                    if ((q != p) && (! gem_node_is_marked(q, nvis, visP->e)))
                      { /* Stack {q}: */
                        gem_ref_vec_expand(visP, nvis);
                        visP->e[nvis] = q; q->label = nvis;  nvis++;
                      }
                  }
              }
            /* Node {visP->e[tvis]} is completely visited now. */
            tvis++;
          }
      }
    /* Return: */
    assert(nvis <= visP->ne);
    (*nvisP) = nvis;
  }

void gem_traverse(gem_ref_t root, int d, gem_ref_vec_t *visP, int *nvisP)
  { demand((d >= 0) && (d <= gem_DIM_MAX), "invalid dimension");
    gem_domain_traverse(root, d+1, NULL, visP, nvisP);
  }

void gem_domains_enum
  ( gem_ref_t root, 
    int nr, int R[], 
    int ns, int S[],
    gem_ref_vec_t *visP,
    int *nvisP,
    gem_ref_vec_t *repP,
    int *nrepP
  )
  { demand(root != NULL, "null {root}");
    demand(visP != NULL, "null {visP}");
    demand(nvisP != NULL, "null {nvisP}"); 
    int nvis = (*nvisP);  /* Number of nodes seen so far (completely or partially). */
    demand((nvis >= 0) && (nvis <= visP->ne), "invalid {*nvisP}");
   
    /* Representatives of {R}-domains: */
    demand(repP != NULL, "null {repP}");
    demand(nrepP != NULL, "null {nrepP}"); 
    int nrep = (*nrepP);  /* Number of representatives already collected. */
    demand((nrep >= 0) && (nrep <= repP->ne), "invalid {*nrepP}");
    
    /* Stack the root vis: */
    if (! gem_node_is_marked(root, nvis, visP->e))
      { int tvis = nvis;/* Number of nodes in {visP} whose {S}-links have been followed. */
        /* Stack the root as a representative, enumerate its {R}-domain: */
        gem_ref_vec_expand(repP, nrep);
        repP->e[nrep] = root; nrep++;
        gem_domain_traverse(root, nr, R, visP, &nvis);
        while (tvis < nvis)
          { /* Grab the next incompletely visited node {p}: */
            gem_ref_t p = visP->e[tvis];;
            /* Check its S-links: */
            int js;
            for (js = 0; js < ns; js++) 
              { int i = (S == NULL ? js : S[js]); 
                if (i <= p->dim)
                  { /* Get the {i}-neighbor {q} of {p}, collect its R-domain if not marked: */
                    gem_ref_t q = p->adj[i];
                    if ((q != p) && (! gem_node_is_marked(q, nvis, visP->e)))
                      { /* Another {R}-domain found: */
                        gem_ref_vec_expand(repP, nrep);
                        repP->e[nrep] = q; nrep++;
                        /* Stack the entire {R}-domain of {q}: */
                        gem_domain_traverse(q, nr, R, visP, &nvis);
                      }
                  }
              }
            /* Node {visP->e[tvis]} is completely visited now. */
            tvis++;
          }
      }
    (*nvisP) = nvis;
    (*nrepP) = nrep;
  }

vec_typeimpl(gem_ref_vec_t,gem_ref_vec,gem_ref_t);

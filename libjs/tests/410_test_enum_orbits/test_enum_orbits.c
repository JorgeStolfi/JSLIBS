#define PROG_NAME "test_enum_orbits"
#define PROG_DESC "tests the {enum_orbits.h} procedures"
#define PROG_VERS "1.1"

/* Last edited on 2024-11-22 21:03:18 by stolfi */
/* Created on 2007-01-31 by J. Stolfi, UNICAMP */

#define PROG_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <ref.h>
#include <affirm.h>
#include <jsmath.h>
#include <ix.h>

#include <enum_orbits.h>

int32_t main (int32_t argc, char **argv);

void test_enum_items_cycle(void);
  /* Tests the {enum_items} and {enum_cycle} procedures. */

void test_enum_orbits(void);
  /* Tests the {enum_orbits} procedure. */
  
/* !!! Should test the other procedures, e.g. {enum_orbits}. !!! */
  
int32_t main (int32_t argc, char **argv)
  { test_enum_items_cycle();
    test_enum_orbits();
    return 0;
  }
  
void test_enum_items_cycle(void)
  { 
    /* The enumeration items will be addresses into the array {obj}: */
    uint32_t nobj;      /* Number of objects in {obj}. */
    uint32_t *obj;   /* Valid objects are {obj[0..nobj-1]}. */
    uint32_t ps;      /* Position of special elment in {obj}. */

    /* Parameters of the regular object array in Fortran order: */
    uint32_t d = 5;                       /* Dimension of array. */
    ix_size_t sz[5] = { 4, 2, 3, 2, 2 }; /* Array size along each axis. */
    ix_pos_t bp = 0;                     /* Position of array elem {(0,.. 0)}. */
    ix_step_t st[d];                     /* Position increments along each axis. */
    ix_order_t ixor = ix_order_F;        /* Order of index variation. */
    ix_packed_steps((ix_dim_t)d, sz, ixor, st);
    
    auto void make_objects(void);
      /* Stores in {obj} the address of a new vector whose elements
        are the objects of the enumeration; and in {nobj} the number
        of those objects, and in {ps} the position of the special
        object.
        
        The objects comprise a regular array, with {sz[i]} elements
        along each axis {i}, linearized in the C/Pascal fashion;
        followed by one /special/ element.  Each element of the array
        contains the corresponding index tuple, encoded as the digits
        of a decimal number, with index {ix[d-1]} being the units
        digit. The special element contains an all-nines integer with
        {d} digits. */
        
    auto ix_pos_t obj_pos(ref_t r);
      /* Given the address of an object, returns its position 
        in the linear vector {obj} (an integer in {0..nobj-1}). */
        
    auto void print_obj(ref_t r);
      /* Prints the object whose address is {r}.  Assumes that {r}
        is the address of an element of {obj[0..nobj-1]} and 
        that it scontains a recognizable {d}-digit decimal integer. */

    /* Fill {obj[0..nobj-1]} with recognizable stuff: */
    make_objects();
    
    /* Define the step and visit procedure tables to use in the tests: */
    uint32_t nt = 4;  /* Number of step functions. */
    enum_step_t *stp[4];
    enum_visit_t *vis[5];
    
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "testing with %d step functions\n", nt);
    
    uint32_t nv; /* The number of visited nodes. */
    uint32_t maxv; /* The max number of visited nodes. */
        
    auto ref_t step0(ref_t r); 
    auto ref_t step1(ref_t r);
    auto ref_t step2(ref_t r);
    auto ref_t step3(ref_t r);
    
    auto ref_t gen_step(ref_t r, ix_axis_t i, uint32_t x);
      /* The general step function: cyclic shift by {x} along axis {i}. */
    
    auto bool_t visit0(ref_t r);
    auto bool_t visit1(ref_t r);
    auto bool_t visit2(ref_t r);
    auto bool_t visit3(ref_t r);
    auto bool_t visit4(ref_t r);
    
    auto bool_t gen_visit(ref_t r, ix_axis_t i);
      /* The generic visit function with index {i}: prints
        "visit{i}({*r})", increments {nv}, and and returns TRUE iff
        {nv >= maxv}. */
   
    assert(nt == 4);
    vis[0] = &visit0; stp[0] = &step0;
    vis[1] = &visit1; stp[1] = &step1;
    vis[2] = &visit2; stp[2] = &step2;
    vis[3] = &visit3; stp[3] = &step3;
    vis[4] = &visit4;
    
    /* The root set is {r[0..2]}: */
    ref_t r[3]; 
    r[0] = &(obj[ 3]); 
    r[1] = &(obj[11]); 
    r[2] = &(obj[79]);
    
    ref_vec_t v = ref_vec_new(0); /* The list of visited items. */
        
    auto void init_visit_list(ref_vec_t *vP);
    /* Initializes the visit list {vP} to a single ref to the letter '*' */
    
    auto void print_visit_list(ref_vec_t *vP);
    /* Prints the visit list {vP} to stderr. */

    for (uint32_t trunc = 0;  trunc < 2; trunc++)
      { /* {trunc = 0} full enum, {= 1} truncated enum. */
        maxv = (trunc == 0 ? nobj+1 : 5);
        
        /* Test {enum_cycle}. */
        for (uint32_t i = 0;  i < nt; i++)
          { uint32_t *r = &(obj[7]);
            fprintf(stderr, "testing enum_cycle with");
            fprintf(stderr, "  step = stp[%d]  visit = vis[%d]", i, i);
            fprintf(stderr, "  root = "); print_obj(r);
            fprintf(stderr, "  maxv = %d", maxv);
            fprintf(stderr, "\n");
            init_visit_list(&v);
            nv = 0;
            enum_cycle((ref_t)r, stp[i], vis[i], &v);
            print_visit_list(&v);
          }

        fprintf(stderr, "testing enum_items");        
        fprintf(stderr, " with stp[0..%d], vis[0..%d]", nt-1, nt-1);        
        fprintf(stderr, "\n");
        init_visit_list(&v);
        nv = 0;
        enum_items(ref_vec_make_desc(r,3), nt, stp, vis, &v);
        print_visit_list(&v);
      }
      
    fprintf(stderr, "============================================================\n");
    return;
    
    /* IMPLEMENTATIONS OF INTERNAL PROCS */
    
    ix_pos_t obj_pos(ref_t r) 
      { uint32_t *ur = (uint32_t *)r;
        assert(obj <= ur);
        uint32_t p = (uint32_t)(ur - obj);
        assert(p < nobj);
        return (uint32_t)p;
      }
      
    /* Step procedures (permutations of {{ &(obj[i]) : i in 0..nobj-1 }}: */
   
    ref_t step0(ref_t r)
      { return gen_step(r, 0, 1); }

    ref_t step1(ref_t r)
      { return gen_step(r, 1, 1); }

    ref_t step2(ref_t r)
      { return gen_step(r, 2, 1); }

    ref_t step3(ref_t r)
      { return gen_step(r, 0, 2); }

    ref_t gen_step(ref_t r, ix_axis_t i, uint32_t x)
      { /* Cycles by {+d} each group of {m} items: */
        ix_pos_t p = obj_pos(r);
        assert(p < nobj);
        /* The object must be part of the regular array: */
        assert(p != ps);
        ix_index_t ix[d];
        ix_packed_indices((ix_dim_t)d, p, bp, sz, ixor, ix);
        assert(ix[i] >= 0);
        ix[i] = (ix_index_t)(((uint32_t)(ix[i]) + x) % sz[i]); 
        p = ix_packed_position((ix_dim_t)d, ix, bp, sz, ixor);
        assert(p < nobj);
        assert(p != ps);
        return (ref_t)(obj + p);
      }

    /* Visit procedures (print the {obj} entry addressed by {r}): */
   
    bool_t visit0(ref_t r)
      { return gen_visit(r, 0); }
      
    bool_t visit1(ref_t r)
      { return gen_visit(r, 1); }

    bool_t visit2(ref_t r)
      { return gen_visit(r, 2); }

    bool_t visit3(ref_t r)
      { return gen_visit(r, 3); }

    bool_t visit4(ref_t r)
      { return gen_visit(r, 4); } 
    
    bool_t gen_visit(ref_t r, ix_axis_t i)
      { fprintf(stderr, "  visit%u(", i);  
        print_obj(r);
        fprintf(stderr, ")\n");
        nv++;
        return nv >= maxv;
      }
    
    void init_visit_list(ref_vec_t *vP)
      { /* Discard all objects, store only the special object: */
        ref_vec_trim(vP, 1);
        vP->e[0] = &(obj[ps]); 
      }

    void print_visit_list(ref_vec_t *vP)
      { fprintf(stderr, "v = [");
        for (uint32_t j = 0;  j < v.ne; j++) { fputc(' ', stderr); print_obj(v.e[j]); }
        fprintf(stderr, " ]\n");
      }

    void print_obj(ref_t r)
      { uint32_t *ru = (uint32_t *)r;
        assert(ru - obj >= 0);
        assert(ru - obj < nobj);
        fprintf(stderr, "%0*u", d, *ru);
      }

    void make_objects(void)
      {
        uint32_t na = (uint32_t)ix_num_tuples((ix_dim_t)d, sz); /* Number of objects in regular array */
        nobj = na + 1;  /* Total number of objects. */
        obj = talloc(nobj, uint32_t);
        /* Fill the regular array elements: */
        ix_index_t ix[d]; 
        assert(ix_assign_min((ix_dim_t)d, ix, sz));
        ix_pos_t p = ix_position((ix_dim_t)d,ix,bp,st);
        uint32_t id;
        uint32_t m = 0; /* Consistency check. */
        do 
          { /* Pack the indices {ix[0..d-1]} as digits of a decimal integer: */
            id = 0;
            for (uint32_t i = 0;  i < d; i++) { assert(ix[i] < 9); id = (uint32_t)(10*id + ix[i]); }
            obj[p] = id;
            /* fprintf(stderr, "p = %llu  m = %d\n", p, m); */
            assert(p == m);
            m++;
          }
        while (! ix_next((ix_dim_t)d,ix,sz,ixor,st,&p,NULL,NULL,NULL,NULL));
        assert(m == na);
        /* Append the special object {obj[ps]}: */
        ps = m; m++;
        id = 0; for (uint32_t i = 0;  i < d; i++) { id = 10*id + 9; }
        obj[ps] = id;
        assert(m == nobj);
      }

  }

void test_enum_orbits(void)
  { 
    /* The enumeration items will be addresses into the array {obj}: */
    uint32_t nobj;      /* Number of objects in {obj}. */
    uint32_t *obj;   /* Valid objects are {obj[0..nobj-1]}. */
    uint32_t ps;      /* Position of special element in {obj}. */
    
    uint32_t nc = 60;

    auto void make_objects(void);
      /* Stores in {obj} the address of a new vector whose elements
        are the objects of the enumeration; and in {nobj} the number
        of those objects, and in {ps} the position of the special
        object.
        
        The objects comprise a vector of {obj[0..nc-1]} elements,
        followed by a special element {obj[ps]}; so {nobj == nc+1} and
        {ps == nc}. Each element of the array contains the
        corresponding index as an integer. The special element
        {obj[ps]} contains an all-nines integer. */
        
    auto ix_pos_t obj_pos(ref_t r);
      /* Given the address of an object, returns its position 
        in the linear vector {obj} (an integer in {0..nobj-1}). */
        
    auto void print_obj(ref_t r);
      /* Prints the object whose address is {r}.  Assumes that {r}
        is the address of an element of {obj[0..nobj-1]} and 
        that it scontains a decimal integer. */

    /* Fill {obj[0..nobj-1]} with recognizable stuff: */
    make_objects();
    
    uint32_t nv; /* The number of visited nodes. */
    uint32_t maxv; /* The max number of visited nodes. */
        
    auto ref_t step2(ref_t r); /* Increment by 2 modulo {nobj}. */
    auto ref_t step3(ref_t r); /* Increment by 3 modulo {nobj}. */
    auto ref_t step4(ref_t r); /* Increment by 4 modulo {nobj}. */
    auto ref_t step5(ref_t r); /* Increment by 5 modulo {nobj}. */
//      auto ref_t step6(ref_t r); /* Increment by 6 modulo {nobj}. */
   
    auto ref_t gen_step(ref_t r, uint32_t x);
      /* The general step function: Increment by {x} modulo {nobj}. */
    
    auto bool_t visit(ref_t r);
      /* The generic visit function with index {i}: prints
        "visit({*r})", increments {nv}, and and returns TRUE iff
        {nv >= maxv}. */
   
    /* Define the {istep}, {ostep} procedure tables to use in the tests: */
    uint32_t ni = 2;  /* Number of intra-subgroup step functions. */
    enum_step_t *istep[2];
    istep[0] = &step5;
    istep[1] = &step4;
    assert(ni == 2);
    
    uint32_t no = 2;  /* Number of inter-subgroup step functions. */
    enum_step_t *ostep[2];
    ostep[0] = &step2;
    ostep[1] = &step3;
    
    /* The root set is {r[0..nr]}: */
    uint32_t nr = 2;
    ref_t r[nr]; 
    r[0] = &(obj[ 0]); 
    r[1] = &(obj[11]); 
    assert(nr == 2);
    
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "testing with %d istep and %d ostep functions\n", ni, no);
    
    ref_vec_t v = ref_vec_new(0); /* The list of visited items. */
        
    auto void init_visit_list(ref_vec_t *vP);
    /* Initializes the visit list {vP} to a single ref to the letter '*' */
    
    auto void print_visit_list(ref_vec_t *vP);
    /* Prints the visit list {vP} to stderr. */

    for (uint32_t trunc = 0;  trunc < 2; trunc++)
      { /* {trunc = 0} full enum, {= 1} truncated enum. */
        maxv = (trunc == 0 ? nobj+1 : 5);
        fprintf(stderr, "testing enum_orbits");        
        fprintf(stderr, " with istep[0..%d], ostep[0..%d]", ni-1, no-1);        
        fprintf(stderr, "\n");
        init_visit_list(&v);
        nv = 0;
        enum_orbits(ref_vec_make_desc(r,nr), ni, istep, no, ostep, &visit, &v);
        print_visit_list(&v);
      }
      
    fprintf(stderr, "============================================================\n");
    return;
    
    /* IMPLEMENTATIONS OF INTERNAL PROCS */
    
    ix_pos_t obj_pos(ref_t r) 
      { uint32_t *ur = (uint32_t *)r;
        assert(obj <= ur);
        uint32_t p = (uint32_t)(ur - obj);
        assert(p < nobj);
        return (ix_pos_t)p;
      }
      
    /* Step procedures (permutations of {{ &(obj[i]) : i in 0..nobj-1 }}: */
   
    ref_t step2(ref_t r)
      { return gen_step(r, 2); }

    ref_t step3(ref_t r)
      { return gen_step(r, 3); }

    ref_t step4(ref_t r)
      { return gen_step(r, 4); }

    ref_t step5(ref_t r)
      { return gen_step(r, 5); }

//      ref_t step6(ref_t r)
//        { return gen_step(r, 6); }

    ref_t gen_step(ref_t r, uint32_t x)
      { /* Cycles by {+d} each group of {m} items: */
        ix_pos_t p = obj_pos(r);
        assert(p < nobj);
        /* The object must be part of the regular array: */
        assert(p != ps);
        p = (p + x) % nc;
        assert(p < nobj);
        assert(p != ps);
        return (ref_t)(obj + p);
      }

    /* Visit procedure (print the {obj} entry addressed by {r}): */
    
    bool_t visit(ref_t r)
      { fprintf(stderr, "  visit(");  
        print_obj(r);
        fprintf(stderr, ")\n");
        nv++;
        return nv >= maxv;
      }
    
    void init_visit_list(ref_vec_t *vP)
      { /* Discard all objects, store only the special object: */
        ref_vec_trim(vP, 1);
        vP->e[0] = &(obj[ps]); 
      }

    void print_visit_list(ref_vec_t *vP)
      { fprintf(stderr, "v = [");
        int32_t j;
        for (j = 0; j < v.ne; j++) { fputc(' ', stderr); print_obj(v.e[j]); }
        fprintf(stderr, " ]\n");
      }

    void print_obj(ref_t r)
      { uint32_t *ru = (uint32_t *)r;
        assert(ru - obj >= 0);
        assert(ru - obj < nobj);
        fprintf(stderr, "%03u", *ru);
      }

    void make_objects(void)
      {
        nc = 60;
        nobj = nc + 1;  /* Total number of objects. */
        obj = talloc(nobj, uint32_t);
        /* Fill the regular array elements: */
        for (uint32_t i = 0; i < nc; i++) { obj[i] = i; }
        /* Append the special object {obj[ps]}: */
        ps = nc;
        obj[ps] = 99;
      }

  }
 

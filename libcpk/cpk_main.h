#ifndef cpk_main_H
#define cpk_main_H

/* Specifying and solving a RadCom problem instance. */
/* Last edited on 2024-12-31 16:16:24 by stolfi */

#include <stdint.h>

#include <vec.h>
#include <r2.h>
#include <bool.h>

#include <cpk_weight.h>

/* 
  THE COMMUNITARY RADIO ALLOCATION PROBLEM
  
  This module is concerned with the problem of allocating community
  radio (CR) stations over a specified region, typically a
  municipality.
  
  All CR stations are assumed to use the same frequency, so the
  distance between any two authorized stations must be at least {dMin
  = 2*rRad}, where {rRad} is the (nominal) maximum broadcast range of
  a CR station. The regulating agency faces the problem ofensuring
  this minimum separation while trying to satisfy the demand for such
  radios by potential operators and potential listeners.
  
  To this end, the agency plans to hold a set of public /auctions/.
  Each auction is held for a specific /auction area/. Among all
  parties interested in installing a CR station within that area, only
  one will be selected (by unspecified criteria) and authorized to
  operate.
  
  Since the auctions are held in parallel, and the auction winner may
  install his station anywhere within the auction area, one must
  observe a distance of at least {dMin = 2*rRad} between any two auction
  areas. It has already been decided that all auction areas will be
  disks of a fixed radius {rAuc}; thus the distance between the
  centers of two such disks must be at least {dMin + 2*rAuc}. Typical
  values are {dMin = 4 km}, {rAuc = 1 km}, which gives a minimum
  distance {4 km} between two auction disks, i.e. {6 km} between two
  auction centers.
  
  There are other geographical constraints to consider. It would be
  silly to open an auction in an unpopulated or inaccessible region,
  such as a lake or nature preserve; thus the center of every auction
  disk must be within the urbanized area. Also, for administrative
  reasons, the auction area must lie completely within a specified
  political region, such as a municipality; and it may be necessary to
  require a minimum distance between any auction area and adjacent
  municipalities. One must also have to observe a (larger) minimum
  distance to other nations.
  
  The problem addressed by this module is therefore to propose an
  optimal set of auctions that honors those constraints. The problem's
  input data is the geography of the region, including the outline of
  the urbanized area, the political boundaries (with other
  municipalities and other nations) that must be honored, existing and
  planned CR stations, and the known demand for such radios.
  
  The known demand is a list of formal `declarations of interest'
  (DIs); each item is a request by some entity for permission to
  install a CR station at a specific location. It should be noted
  that, once the auctions are officially opened, any party may enter
  the contest for a given area, independently of whether they have
  previously submitted a DI or not. Thus the DIs should be viewed
  only as rough estimates of the actual demand.

  The output of the problem is a set of proposed auction disks, which
  obey the minimum separation constraint and are as good as possible
  under various criteria --- such as how many auctions would be
  opened, how many of the existing DIs would be eligible to enter
  them, and how much of the urbanized area would be covered by the
  winning stations.
  
  In this package, the auction centers are drawn from some finite list
  of candidate points {V[0..nV]} of {\Real^2}, typically belonging to
  some regular discrete grid. Thus a potential solution will be
  internally represented as a list {J[0..nJ-1]} of indices into {V},
  in the range {0..nV-1}. */
  
/*
  GEOGRAPHIC DATA RECORD
  
  The input data is packaged as a single record of type {cpk_domain_t},
  defined below. */

typedef struct cpk_domain_t  /* Specifies the valid points of {\Real^2} */
  { r2_vec_t Urb;       /* Vertices of urbanized area polygon. */
    r2_vec_t Exs;       /* Existing or planned stations. */
    r2_vec_t Auc;       /* Centers of predetermined auction zones. */
    r2_vec_t Mun;       /* Vertices of intermunicipal borders. */
    r2_vec_t Nat;       /* Vertices of international borders. */
    r2_vec_t Dem;       /* Known demand points (DI)s. */
  } cpk_domain_t;  
  /* 
    A {cpk_domain_t} {C} specifies a subset of the plane {\Real^2},
    the /{C}-valid/ points, that can be considered as auction centers;
    and the declared demand in that area.
    
    Given the policy parameters {dMin,rAuc,dMun,dNat} (see below), a
    point {x} is {C}-valid iff
    
      (0) {x} lies inside the polygon {Urb}; and
      (1) {dist(x,s) >= dMin + rAuc} for any station points {s} in {Exs}; and
      (2) {dist(x,a) >= dMin + 2*rAuc} for any auction center {a} in {Auc}; and
      (3) {dist(x,m) >= dMun} for any point {m} on the polyline {Mun}; and
      (4) {dist(x,n) >= dNat} for any point {n} on the polyline {Nat}.
    
    Typically, {Urb} would be the urbanized area, {Mun} the boundaries
    of adjacent municipalities (or, more generally, of any areas where
    CR stations could be independently allocated), and {Nat} the
    boundaries of other nations. Each point in {Exs} would be the
    location of an existing station, and each point in {Auc} would be
    the center of an auction area already decided upon.
    
    The polygon {Urb} may have holes an multiple components; see
    {hxg_canvas_polygon} for details of its representation. Similarly
    the the polylines {Mun} and {Nat} may have multiple separated
    pieces; see {hxg_canvas_polyline} for details.
    
    All of the vectors {Exs,Auc,Mun,Nat,Dem} may be empty, but {Urb}
    must always be provided. (If the real urbanized area is not known,
    the client must provide some surrogate closed polygon, e.g. the
    municipality outline.) */

void cpk_domain_free(cpk_domain_t *C);
  /* Calls {free} on all vectors hanging from {C}. */

/* 
  POLICY PARAMETERS 
  
  In addition to the geographic data, the problem's input includes a
  set of /policy parameters/ defined by the final user (such as the
  regulating agency or local authorities). These include the radius
  {rAuc} of the auction areas, various minimum distances to be
  observed, and the formula to be used when comparing potential
  solutions. These parameters are packaged as a {cpk_policy_t}
  record, defined below. */

typedef struct cpk_policy_t
  { double dMin;    /* Minimum distance between radio stations. */
    double rAuc;    /* Radius of auction areas. */
    double dMun;    /* Minimum distance from auction center to other municipalities. */
    double dNat;    /* Minimum distance from auction center to other nations. */
    bool_t for_demand;     /* Auctions must be near DIs, or fill whole urban area? */
    cpk_wcoeffs_t wcoeffs; /* Coefficients of the weight formula. */
  } cpk_policy_t;
  /* Policy parameters defining the constraints and ranking of 
    solutions. 

    If {for-demand} is TRUE, every proposed auction disk will include
    at least one DI. Said another way, the auction centers are
    restricted to the area within distance {rAuc} of some {DI}. In
    particular, if there are no DIs, the output will be empty.

    If {for-demand} is FALSE, the proposed solution will typically
    fill the whole urban area {C->Urb}, even if there are no DIs.
    The DIs may still affect the ranking of such solutions. */

void cpk_policy_free(cpk_policy_t *P);
  /* Calls {free} on all vectors hanging from {P}. */
  
/*
  OPTIMIZATION OPTIONS
  
  Finally, the computation of a good auction set is controlled by
  a set of /options/ which can be set by the user and/or the 
  calling program. They do not affect the set of valid solutions
  nor the criteria for comparing them, but only /how/ the procedure
  will look for good solutions, and how much time it can spend on 
  the search. They also control the amount of diagnostic messages
  and plot files, and internal consitency checks. These options
  are packaged as a {cpk_options_t} record defined below.  */

typedef struct cpk_options_t
  { bool_t verbose;   /* TRUE prints various diagnostics. */
    bool_t validate;  /* TRUE runs (expensive) validation checks on solution. */
    bool_t plot;      /* TRUE plots the main solutions (greedy and grasp). */
    uint32_t maxSecsGRASP; /* Maximum number of seconds allowed for GRASP. */
    double magnify;   /* Domain enlargement factor (for testing). */
    uint32_t seed;         /* A seed for random choices. */
    char *outDir;     /* Directory for output file names. */
  } cpk_options_t;

/* 
   MAIN PROCEDURE */

void cpk_choose_auction_points 
  ( cpk_domain_t *C,    /* Geographic data, exs/pln stations, and declared demand ("DI"s). */
    cpk_policy_t *P,    /* Allocation policy parameters. */
    cpk_options_t *opt, /* Algorithm options. */
    r2_vec_t *VJ,       /* (OUT) Proposed solution (list of auction points). */
    double *WJ          /* (OUT) Score of proposed solution. */
  );
  /* 
    Returns in {VJ} a list of proposed auction centers ("pontos
    simulados") for the region described by {C}, according to the
    parameters in {P}. Also returns in {WJ} the numerical score of
    that list, which it tries to maximize.
    
    All coordinates (input and output) are Longitude/Latitude pairs,
    in fractional degrees. They are converted internally to/from EUTM
    X/Y variables (in meters), as described in the documentation of
    {cpk_LL_to_EUTM} and {cpk_EUTM_to_LL}; and then multiplied by
    {opt->magnify} (usually 1.0). The reference longitude {refLon} and
    reference EUTM-Y {refY} used in the conversion are selected by the
    procedure (approximate mean coordinates of the urbanized area). */

#endif

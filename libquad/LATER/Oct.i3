(* Oct-edge data structure. *)
(* Created 1993 by Rober Marcone Rosi and J. Stolfi *)
(* See 

   "Primitives for the Manipulation of General Subdivisions 
   and the Computation of Voronoi Diagrams"

   L. Guibas, J. Stolfi, ACM TOG, April 1985
*)

INTERFACE Oct;

(* A 'octet data structure' is a pair of dual subdvisions of a connected
   2-D manifold.
   
   An 'edge' is an edge of the primal or dual subdivision.
   
   A 'dired' is a directed edge.
   
   An 'ored' is an oriented edge, i.e. an edge with local orientation
   (a 'side') of the manifold.
   
   A 'quad' is a group of four arcs, corresponding to two mutually dual
   edges, each taken in its two possible senses and with a fixed local
   orientation.
   
   An 'octet' is a group of eight arcs, consisting of a pair of mutually 
   dual edges in all their directions and local orientations.
*)

IMPORT Wr, Rd;

TYPE
  T = Arc; (* A subdivision is usually handled by one of its arcs. *)

  Arc = RECORD edge: Edge; bits: FRBits END; (* A dir'd+ori'd edge *)
  
  FRBits = [0..7]; (* Direction and side w.r.t. base arc. *)
  SBit = [0..1];   (* Sense bit *)
  FBit = [0..1];   (* Flip bit *)
  RBits = [0..3];  (* Rotation bits *)
  DBit = [0..1];   (* Dual/primal bit *)
  
  ArcNum = [0..MaxArcNum];
  EdgeNum = [0..MaxEdgeNum];

  Edge <: PublicEdge;
  PublicEdge = OBJECT 
      num: EdgeNum := 0;
    METHODS
      init(): Edge;
        (* Initializes the pointers of "self" so that it
           describes an isolated edge with distinct endpoints.
           That means a subdivision of the sphere with one edge,
           two vertices, and one face. Returns "self". *)
    END;
  
(* ====== Arcs ====== *)

(* A "FRBits" value identifies a particular arc within any edge, according 
   to the following table (where 'e' is the reference arc of the edge):
    
|     FRBits    arc            FBit  SBit  DBit   RBits
|     --------------------------------------------------
|      000 0    e               0     0     0       0  
|      001 1    e.flip          1     0     0       0  
|      010 2    e.rot           0     0     1       1  
|      011 3    e.rot.flip      1     0     1       1  
|      100 4    e.sym           0     1     0       2  
|      101 5    e.sym.flip      1     1     0       2  
|      110 6    e.tor           0     1     1       3  
|      111 7    e.tor.flip      1     1     1       3   
      
   In general, the arc "Flip^f(Rot^r(e))", for "r" in [0..3] and "f" in [0..1],
   has "FRBits" equal to "2r + f".
*)

(* ====== Edge orientation operators ====== *)

PROCEDURE Rot  (a: Arc): Arc; (* dual edge directed from right to left face *)
PROCEDURE Flip (a: Arc): Arc; (* reverses left and right faces *)
PROCEDURE Sym  (a: Arc): Arc; (* "a Sym = a Rot Rot" *)
PROCEDURE Tor  (a: Arc): Arc; (* Inverse of "Rot" *)

(* ====== Vertex/face walking operators ====== *)

PROCEDURE Onext (a: Arc): Arc; (* next arc in ccw order with same origin *)
PROCEDURE Dnext (a: Arc): Arc; (* next arc in ccw order with same dest. *)
PROCEDURE Rnext (a: Arc): Arc; (* next arc in ccw order with same rgt face *) 
PROCEDURE Lnext (a: Arc): Arc; (* next arc in ccw order with same lft face *)

PROCEDURE Oprev (a: Arc): Arc; (* "a Onext^(-1)" *) 
PROCEDURE Dprev (a: Arc): Arc; (* "a Dnext^(-1)" *)
PROCEDURE Rprev (a: Arc): Arc; (* "a Rnext^(-1)" *)
PROCEDURE Lprev (a: Arc): Arc; (* "a Lnext^(-1)" *) 

PROCEDURE Degree(a: Arc): CARDINAL;
  (* Number of "Arc"s with same origin as "a" *)

PROCEDURE SenseBit (a: Arc): SBit;  
  (* Sense bit = "FRBits DIV 4".  
     Reversed by "Sym", unchanged by "Flip", may be changed by "Rot". *)

PROCEDURE FlipBit (a: Arc): FBit;   
  (* Flip bit = "FRBits MOD 2".
     Reversed by "Flip", unchanged by "Sym" and "Rot". *)

PROCEDURE DualBit (a: Arc): DBit;
  (* Dual bit = "FRBits DIV 2 MOD 2".
     Unchanged by "Flip" and "Sym", reversed by "Rot". *)

PROCEDURE RotBits (a: Arc): RBits;
  (* Rot bits = "FRBits DIV 2".
     Unchanged by "Flip", incremented twice by "Sym",
     either incremented or decremented by "Rot" (depending
     on the "Flip" bit). *)

(* ====== Creation ====== *)

PROCEDURE MakeEdge (): Arc;
  (* Creates a new edge. Equivalent to 
  |     Arc{edge := NEW(Edge).init(), bits := 0}
  *)

(* ====== Splicing ====== *)

PROCEDURE Splice (a, b: Arc);
(* Merges or splits the origin and left faces of a and b.
   Assumes "a" cannot be reached from "b" through an
   odd number of "Rot"'s, and also "a" is not "Flip(Onext(b))". *)
     
PROCEDURE SetOnext(a, b: Arc);
(* If "Onext(a) # b", performs "Splice(a, Oprev(b))".  
   After this call, "Onext(a)" will be equal to "b".
   Valid whenever "Splice(a,b)" is valid. *)

(*
  Conjecture: let "e[0..n]" be all arcs of an "Oct" structure
  with "FlipBit = 0"; and "x[0..n]" a list of arcs 
  such that "Onext(e[i]) = x[i]".
  
  We can build a copy of the structure by the following
  algorithm:
  
    1. create n/4 new "Edge" records, and a vector "ec[0..n]" of
       their unFlipped arcs, such that
         (a) if "Rot(e[i]) = e[j]", then "Rot(ec[i]) = ec[j]"
         
    2. Create a vector "xc[0..n]" of Arc's, such that if 
       "x[i] = Flip^f(e[j])", then "xc[i] = Flip^f(ec[j])".
    
    3. for "i = 0..n", in any order, perform "SetOnext(ec[i], xc[i])".
    
  This conjecture is true if we can prove that 
  
    every arc has its "Onext" properly set;
    
    whenever "SetOnext(ec[i], ec[x[i]])" changes "Onext(e[j])"
    for some arc "j", then we still haven't performed 
    "SetOnext(ec[j], ec[x[j]])".
    
*)

(* ====== Traversal ====== *)

TYPE 
  VisitProc = PROCEDURE (a:Arc);
    (* A client-provided arc visitation procedure *)
  
PROCEDURE EnumEdges(a: Arc; visit: VisitProc; edges: BOOLEAN := FALSE); 
  (* Enumerates all arcs  that can be reached from a by some combination 
     of Sym's and Onext's.

     If "edges" is "FALSE" (default), each arc is passed to the VisitProc 
     exactly once.

     If "edges" is "TRUE", each UNDIRECTED edge is passed to the VisitProc
     exactly once, in one of its two directions.

     If the manifold is orientable, then each arc will be visited in 
     only one of its orientations. Otherwise, for every visited arc "a" the 
     procedure will also visit "Flip(a)" and/or "Flip(Sym(a))".

     If the algebra is well-formed, only arcs of same primduality as a 
     will ever be enumerated; that is, the procedure will never visit 
     "Rot(a)" or "Flip(Rot(a))" if it visits the arc "a". 
  *)
     
PROCEDURE EnumVertices(a: Arc; visit: VisitProc);
  (* Enumerates all (primal) vertices reachable from "a" by chains
     of "Sym" and "Onext". The "VisitProc" will be called on exactly
     one arc out of each vertex. *)
     
(* ====== Numbering ====== *)

CONST 
  MaxArcNum = LAST(CARDINAL);    (* Maximum arc number *)
  MaxEdgeNum = MaxArcNum DIV 8;  (* Maximum edge number *)

PROCEDURE NumberEdges(READONLY a: ARRAY OF Arc): REF ARRAY OF Arc;
  (* Assigns distinct numbers to all edges reachable
     from "a" by Sym/Onext chains. Returns a vector 
     with one reachable arc from each edge. *)

PROCEDURE GetArcNum (a: Arc): ArcNum;
  (* = "a.edge.num * 8 + a.bits" *)

(* ====== Printout ====== *)

PROCEDURE PrintArc(wr: Wr.T; a: Arc; eWidth: CARDINAL := 1);
  (* Prints arc "a" as "n:r:f", where "n" is the serial edge number, and "r,f" 
     are such that "a = Flip^f(Rot^r(a0))", where "a0" is the base arc of the
     edge.  The edge number will be left-padded with spaces to "eWidth" bytes. *)

PROCEDURE ReadArc(rd: Rd.T; READONLY map: ARRAY OF Edge): Arc;
  (* Reads from "rd" an arc in the "n:r:f" format used by "PrintArc".
     The "map" table is used to convert the edge number "n" into 
     an "Edge" pointer. *)

PROCEDURE PrintEdge(wr: Wr.T; e: Edge; eWidth: CARDINAL := 1);
  (* Prints on "wr" the four arcs Onext(Rot^i(s)),
     where "s = Arc{e, 0}" and "i = 0..3".  
     Arcs are separated by one space.
     Each arc is printed using "PrintArc", with edge numbers
     left-padded to "eWidth" bytes. *) 

PROCEDURE ReadEdge(rd: Rd.T; e: Edge; READONLY map: ARRAY OF Edge);
  (* Reads from "rd" four arcs "a[0..3]", using "ReadArc(rd, map)".
     Then performs "SetOnext(Rot^i(s), a[i])",
     where "s = Arc{0,e}" and "i = 0..3". *)

END Oct.













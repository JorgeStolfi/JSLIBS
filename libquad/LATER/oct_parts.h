#ifndef oct_parts_H
#define oct_parts_H

/* Procedures that build oct-edge structures for various simple maps. */
/* Last edited on 2007-01-16 19:12:46 by stolfi */

#define oct_parts_H_copyright \
  "Copyright © 1996, 2006 Institute of Computing, Unicamp."

// MODULE MakeShape EXPORTS Main;
// 
// IMPORT 
//   Scan, Fmt, ParseParams, Text, Oct, Map,
//   Triang, Color, Random, Wr, Stdio, Thread, RTMisc;
// FROM Oct IMPORT 
//   Splice, Enum, Onext, Oprev, Lnext, Rnext, 
//   Rprev, Rot, Sym, Tor, Flip, PrintArc;
// 
// FROM Map IMPORT Arc, GluePatch, Middle;
// 
// FROM Triang IMPORT Org, SetNodeColor, SetVertexRadius;
// 
// FROM Stdio IMPORT stderr;
// 
// TYPE 
//   Shape = {
//     Torus, BiTorus, TriTorus, Klein, Klein2, Klein3, PPlane, 
//     Tetra, Stick, Ring, Cube, Sausage, Orange, Fork, Star
//   };
//   Options = RECORD
//       gridOrder: CARDINAL;
//       shape: Shape;
//       shapeName: TEXT;
//       colr1, colr2: Color.T;
//     END;
// 
// PROCEDURE putwr (a: Arc) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     Wr.PutText(stderr, "edge: "); PrintArc(stderr, a);
//     Wr.PutText(stderr, "  onext: "); PrintArc(stderr, Onext(a));
//     Wr.PutText(stderr, "\n");
//   END putwr;
// 
// PROCEDURE prt (a: Arc) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     putwr(a);
//   END prt;
// 
// PROCEDURE CalcTriang(a: Arc): Triang.Arc =
//   BEGIN
//     Enum(a, GluePatch);
//     Enum(a, prt);
//     RETURN Middle(a)
//   END CalcTriang;
// 
// VAR EdgeCount: CARDINAL := 0;
//   
// (* Basic tools: *)
// 
// PROCEDURE MkEdge(gridOrder: CARDINAL): Arc = 
//   BEGIN
//     WITH 
//       color1 = Color.T{1.00, 0.95, 0.70},
//       color2 = Color.T{1.00, 0.95, 0.70},
//       a = Map.Make(gridOrder, EdgeCount)
//     DO 
//       Map.SetPrimalProperties(
//         a, 
//         vertexRadius := 0.03,
//         vertexColor := Color.T{0.0, 0.0, 0.0},
//         edgeRadius := 0.01,
//         edgeColor := Color.T{0.0, 0.0, 0.0},
//         faceColor := color1,
//         faceTransp := Color.T{0.8, 0.8, 0.8}
//       );
//       Map.SetPrimalProperties(
//         Rot(a),
//         vertexRadius := 0.02,
//         vertexColor := Color.T{1.0, 0.2, 0.0},
//         edgeRadius := 0.01,
//         edgeColor := Color.T{1.0, 0.2, 0.0},
//         faceColor := color2,
//         faceTransp := Color.T{0.8, 0.8, 0.8}
//       );
//       WITH vc = Org(Middle(a)) DO
//         SetNodeColor(vc, Color.T{0.5, 0.1, 0.0});
//         SetVertexRadius(vc, 0.02)
//       END;
//       INC(EdgeCount, 1); 
//       RETURN a
//     END
//   END MkEdge;
// 
// PROCEDURE MakeRing (n: INTEGER; gridOrder: CARDINAL): Arc =
//   VAR fst, a, b: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MkEdge(gridOrder); 
//     fst := a;
//     FOR i := 2 TO n DO
//       b := MkEdge(gridOrder); 
//       Splice(b, Sym(a));
//       a := b;
//     END;
//     Splice(fst, Sym(a));
//     RETURN fst
//   END MakeRing;
// 
// PROCEDURE MakeOrange(n: CARDINAL; gridOrder: CARDINAL): Arc =
//   (* 
//     Same as Rot(MakeRing(n, gridOrder)), except that the 
//     primal/dual coloring of patches is reversed. *)
//   VAR fst, a, b: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := Rot(MkEdge(gridOrder));
//     fst := a;
//     FOR i := 2 TO n DO
//       b := Rot(MkEdge(gridOrder));
//       Splice(b, Sym(a));
//       a := b;
//     END;
//     Splice(fst, Sym(a));
//     RETURN Tor(fst)
//   END MakeOrange;
// 
// (* Shape builders: *)
// 
// PROCEDURE MakeTetra(gridOrder: CARDINAL): Arc =
//   VAR s, t, a: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     s := MakeRing(4, gridOrder);
//     a := MkEdge(gridOrder);
//     t := Rnext(Rnext(s));
//     
//     Splice(s, a); Splice(t, Sym(a));
//     a := MkEdge(gridOrder);
//     Splice(Sym(s), a); Splice(Sym(t), Sym(a));
// 
//     Wr.PutText(stderr, "Exit MakeTetra:\n"); 
//     RETURN a    
//   END MakeTetra;
// 
// PROCEDURE MakeStick (gridOrder: CARDINAL):Arc =
//   VAR a: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MkEdge(gridOrder);
//     Splice(a, Sym(a));
//     Wr.PutText(Stdio.stderr, "Exit MakeStick:\n");
//     RETURN a;
//   END MakeStick;
// 
// PROCEDURE MakeCube (gridOrder: CARDINAL): Arc =
//   VAR b, t, e: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     b := MakeRing(4,gridOrder);
//     t := MakeRing(4,gridOrder); t := Oprev(t);
//     FOR i := 1 TO 4 DO
//       e := MkEdge(gridOrder);
//       Splice(b, e); Splice(Sym(e), t);
//       b := Rprev(b); t := Rnext(t)
//     END;
//     Wr.PutText(stderr, "Exit MakeCube:\n"); 
//     RETURN b    
//   END MakeCube;
// 
// PROCEDURE MakeSausage (length: CARDINAL; gridOrder: CARDINAL): Arc =
//   VAR s, t: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     t := MakeRing(2, gridOrder);
//     BuildTower(2, length, gridOrder, t);
//     Wr.PutText(stderr, "Exit MakeSausage:\n"); 
//     RETURN s    
//   END MakeSausage;
// 
// PROCEDURE MakeFork (prongs, length: CARDINAL; gridOrder: CARDINAL): Arc =
//   VAR o: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     o := MakeOrange(prongs, gridOrder);
//     FOR k := 1 TO prongs DO 
//       o := Onext(o); 
//       BuildTower(2, length, gridOrder, o);
//     END;
//     Wr.PutText(stderr, "Exit MakeFork:\n"); 
//     RETURN o    
//   END MakeFork;
//   
// PROCEDURE BuildTower(m, h, gridOrder: CARDINAL; a: Arc) =
//   (* 
//     Builds a cylindrical tower on the face Right(a), 
//     which must have "m" edges. The tower will have "h" stages
//     and a roof; each stage will be a ring with "m" square faces.
//     The roof of the tower will be a single face of "m" edges. *)
//   VAR s, e, t: Arc;
//   BEGIN
//     t := a;
//     FOR i := 1 TO h DO
//       t := Oprev(t);
//       s := MakeRing(m, gridOrder);
//       FOR j  := 1 TO m DO 
//         e := MkEdge(gridOrder);
//         Splice(t, e); Splice(Sym(e), s);
//         t := Lnext(t); s := Rnext(s)
//       END;
//       t := s;
//     END;
//   END BuildTower;
//   
// PROCEDURE MakePPlane (gridOrder: CARDINAL): Arc =
//   VAR e: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *> 
//   BEGIN
//     e := MkEdge(gridOrder);
//     Splice (Flip(Sym(e)), e);
//     Wr.PutText(stderr, "Exit MakePPlane: \n");
//     RETURN e
//   END MakePPlane;
// 
// PROCEDURE MakeTorus(gridOrder: CARDINAL): Arc =
//   VAR a, b: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MkEdge(gridOrder);
//     b := MkEdge(gridOrder);
//     Splice(a, b);
//     Splice(Sym(a), a);
//     Splice(Sym(b), a);
//     Wr.PutText(stderr, "Exit MakeTorus: \n");    
//     RETURN a
//   END MakeTorus;
//  
// PROCEDURE MakeBiTorus(gridOrder: CARDINAL): Arc =
//   VAR a, b, c, d, e, f: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MkEdge(gridOrder);
//     b := MkEdge(gridOrder);
//     c := MkEdge(gridOrder);
//     d := MkEdge(gridOrder);
//     e := MkEdge(gridOrder);
//     f := MkEdge(gridOrder);
//     Splice(a, b);
//     Splice(b, c);
//     Splice(c, d);
//     Splice(d, e);
//     Splice(e, Sym(c));
//     Splice(Sym(a), Sym(e));
//     Splice(Sym(e), Sym(f));
//     Splice(Sym(f), Sym(d));
//     Splice(Sym(d), Sym(b)); 
//     Splice(Sym(b), f); 
//     Wr.PutText(stderr, "Exit MakeBiTorus: \n");    
//     RETURN a
//   END MakeBiTorus;
//  
// PROCEDURE MakeTriTorus(gridOrder: CARDINAL): Arc =
//   VAR a: ARRAY [0..1] OF ARRAY [0..3] OF Arc;
//       t, s, e: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     (* 
//       Build two tetrahedra, inner and outer.  Store in a[i, k]
//       one arc such that Left(a[i,k]) is face "k" of tetrahedron "i".
//     *)
//     FOR i := 0 TO 1 DO 
//       a[i,0] := MakeTetra(gridOrder);
//       a[i,1] := Onext(a[i,0]);
//       a[i,2] := Sym(Onext(Sym(a[i,0])));
//       a[i,3] := Oprev(a[i,2]);
//     END;
//     (*
//       Build tubes connecting corresponding faces of the two tetrahedra.
//       Each tube has three edges and three square faces.
//     *)
//     FOR k := 0 TO 3 DO 
//       t := a[0,k]; s := a[1,k];
//       FOR j := 1 TO 3 DO 
//         e := MkEdge(gridOrder);
//         Splice(t, e); Splice(s, Flip(Sym(e)));
//         t := Lnext(t); s := Lnext(s)
//       END;
//     END;
//     Wr.PutText(stderr, "Exit MakeTriTorus: \n");    
//     RETURN s
//   END MakeTriTorus;
//  
// PROCEDURE MakeKlein(gridOrder: CARDINAL): Arc =
//   VAR a, b: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MkEdge(gridOrder); 
//     b := MkEdge(gridOrder);
//     Splice(a, b);
//     Splice(Sym(a), a);
//     Splice(Flip(Sym(b)), a);
//     Wr.PutText(stderr, "Exit MakeKlein: \n");
//     RETURN a
//   END MakeKlein;
//   
// PROCEDURE MakeKlein2(gridOrder: CARDINAL): Arc =
//   VAR a, b, c, d: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MkEdge(gridOrder); 
//     b := MkEdge(gridOrder);
//     c := MkEdge(gridOrder); 
//     d := MkEdge(gridOrder);
//     
//     (* Vertex 1: *)
//     Splice(b, a);
//     Splice(a, Flip(Sym(c)));
//     Splice(Flip(Sym(c)), Sym(a));
//     
//     (* Vertex 2: *)
//     Splice(c, d);
//     Splice(d, Sym(b));
//     Splice(Sym(b), Sym(d));
//     
//     Wr.PutText(stderr, "Exit MakeKlein2: \n");
//     RETURN a
//   END MakeKlein2;
//   
// PROCEDURE MakeKlein3(gridOrder: CARDINAL): Arc =
//   VAR a, b, c, d, e, f: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MkEdge(gridOrder); 
//     b := MkEdge(gridOrder);
//     c := MkEdge(gridOrder); 
//     d := MkEdge(gridOrder);
//     e := MkEdge(gridOrder);
//     f := MkEdge(gridOrder);
//     Wr.PutText(stderr, "Exit MakeKlein3: \n");
//     
//     (* Vertex 1: *)
//     Splice(b, a);
//     Splice(a, Flip(Sym(c)));
//     Splice(Flip(Sym(c)), Sym(a));
//     
//     (* Vertex 2: *)
//     Splice(f, d);
//     Splice(d, Sym(b));
//     Splice(Sym(b), Sym(d));
//     
//     (* Vertex 3: *)
//     Splice(c, e);
//     Splice(e, Sym(f));
//     Splice(Sym(f), Sym(e));
//     
//     RETURN a
//   END MakeKlein3;
//   
// PROCEDURE MakeMap(shape: Shape; gridOrder: CARDINAL): Arc =
//   BEGIN
//     CASE shape OF 
//     | Shape.Torus => RETURN MakeTorus(gridOrder);
//     | Shape.BiTorus => RETURN MakeBiTorus(gridOrder);
//     | Shape.TriTorus => RETURN MakeTriTorus(gridOrder);
//     | Shape.Klein => RETURN MakeKlein(gridOrder);
//     | Shape.Klein2 => RETURN MakeKlein2(gridOrder);
//     | Shape.Klein3 => RETURN MakeKlein3(gridOrder);
//     | Shape.PPlane => RETURN MakePPlane(gridOrder);
//     | Shape.Tetra => RETURN MakeTetra(gridOrder);
//     | Shape.Stick => RETURN MakeStick(gridOrder);
//     | Shape.Ring => RETURN MakeRing(5, gridOrder);
//     | Shape.Cube => RETURN MakeCube(gridOrder);
//     | Shape.Sausage => RETURN MakeSausage(3, gridOrder);
//     | Shape.Orange => RETURN MakeFork(7, 0, gridOrder);
//     | Shape.Fork => RETURN MakeFork(3, 2, gridOrder);
//     | Shape.Star => RETURN MakeFork(5, 1, gridOrder);
//     END
//   END MakeMap;
//   
// PROCEDURE Main() =  
//   <* FATAL Wr.Failure, Thread.Alerted *>  
//   BEGIN
//     WITH 
//       o = GetOptions(),
//       m = MakeMap(o.shape, o.gridOrder),
//       s = CalcTriang(m),
//       t = Triang.MakeTopology(s)
//     DO
//       Triang.PrintTopTri(t, o.shapeName & "-" & Fmt.Int(o.gridOrder))
//     END
//   END Main;
// 
// PROCEDURE GetOptions (): Options =
//   <* FATAL Thread.Alerted, Wr.Failure *>
//   VAR o: Options;
//   BEGIN
//     TRY
//       ParseParams.BeginParsing(stderr);                         
// 
//       ParseParams.GetKeyword("-gridOrder");                               
//       o.gridOrder := ParseParams.GetNextInt(1, 20); 
//       
//       ParseParams.GetKeyword("-shape");  
//       o.shapeName := ParseParams.GetNext();
//       IF Text.Equal(o.shapeName, "torus") THEN
//         o.shape := Shape.Torus
//       ELSIF Text.Equal(o.shapeName, "biTorus") THEN
//         o.shape := Shape.BiTorus
//       ELSIF Text.Equal(o.shapeName, "triTorus") THEN
//         o.shape := Shape.TriTorus
//       ELSIF Text.Equal(o.shapeName, "klein") THEN
//         o.shape := Shape.Klein
//       ELSIF Text.Equal(o.shapeName, "klein2") THEN
//         o.shape := Shape.Klein2
//       ELSIF Text.Equal(o.shapeName, "klein3") THEN
//         o.shape := Shape.Klein3
//       ELSIF Text.Equal(o.shapeName, "pPlane") THEN
//         o.shape := Shape.PPlane
//       ELSIF Text.Equal(o.shapeName, "tetra") THEN
//         o.shape := Shape.Tetra
//       ELSIF Text.Equal(o.shapeName, "stick") THEN
//         o.shape := Shape.Stick
//       ELSIF Text.Equal(o.shapeName, "ring") THEN
//         o.shape := Shape.Ring
//       ELSIF Text.Equal(o.shapeName, "cube") THEN
//         o.shape := Shape.Cube
//       ELSIF Text.Equal(o.shapeName, "sausage") THEN
//         o.shape := Shape.Sausage
//       ELSIF Text.Equal(o.shapeName, "orange") THEN
//         o.shape := Shape.Orange
//       ELSIF Text.Equal(o.shapeName, "fork") THEN
//         o.shape := Shape.Fork
//       ELSIF Text.Equal(o.shapeName, "star") THEN
//         o.shape := Shape.Star
//       ELSE
//         Wr.PutText(stderr, "Bad shape \"" & ParseParams.GetNext() & "\"\n");
//         RAISE Scan.BadFormat
//       END;
//       
//       ParseParams.EndParsing();                                       
//     EXCEPT                                                            
//     | Scan.BadFormat =>                                              
//         Wr.PutText(stderr, "Usage: MakeShape -gridOrder <num>\\\n");
//         Wr.PutText(stderr, "  -shape { torus | klein | ... | star }\n");
//         RTMisc.Exit (1);                                              
//     END;
//     RETURN o
//   END GetOptions;
//  
// BEGIN 
//   Main();
// END MakeShape.
//  
// 

// MODULE TestOctet EXPORTS Main;
// 
// IMPORT Oct, SOct, DOct, Wr, Stdio, Thread, Fmt, FileStream;
// 
// FROM Oct IMPORT Splice, Enum, Onext, Oprev, Lnext, Rnext, 
//                 Rprev, Sym, Rot, Flip, PrintArc;
// 
// FROM SOct IMPORT Arc, MakeSOct, GluePatch, Canto, Midle;
// 
// FROM DOct IMPORT Topology, PutRayTriangles, PutWireTriangles, 
//                  MakeTopology, Coords;
// 
// FROM OptShape IMPORT SimplexOpt, GradOpt, PraxOpt;
// 
// FROM GradTest IMPORT GradTOpt, PraxTOpt;
// 
// FROM Stdio IMPORT stderr;
// 
// VAR tetra, ramid, sausage, stick, ring, cube, torus, klein, pplane: Arc;
// 
// CONST GridOrder = 3; (* Order of grid used to realize each quad *)
// 
// PROCEDURE putwr (a: Arc) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     Wr.PutText(stderr, "edge: "); PrintArc(stderr, a);
//     Wr.PutText(stderr, "  onext: "); PrintArc(stderr, Onext(a));
//     Wr.PutText(stderr, "\n");
//   END putwr;
// 
// PROCEDURE prt (a: Arc) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     WITH
//       b = Flip(a), c = Rot(a), d = Flip(c), e = Rot(c),
//       f = Flip(e), g = Rot(e), h = Flip(g)
//     DO
//       putwr(a);
//  (*     putwr(b);
//       putwr(c);
//       putwr(d);
//       putwr(e);
//       putwr(f);
//       putwr(g);
//       putwr(h);  *)
//     END
//   END prt;
// 
// PROCEDURE WriteWireFrame(a: Arc; READONLY coords: Coords; shape: TEXT; tag: TEXT) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     WITH 
//       name = shape & "-" & tag & ".poly",
//       wr = FileStream.OpenWrite(name)
//     DO
//       PutWireTriangles(wr, Rot(a), coords);
//       Wr.Close(wr)
//     END
//   END WriteWireFrame;
// 
// PROCEDURE WriteRayTrace(a: Arc; READONLY coords: Coords; shape: TEXT; tag: TEXT) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     WITH 
//       name = shape & "-" & tag & ".ray",
//       wr = FileStream.OpenWrite(name)
//     DO
//       Wr.PutText(wr, "screen 256 256 \n");       
//       Wr.PutText(wr, "sample 1 \n");       
//       Wr.PutText(wr, "name surfac list \n");       
//       PutRayTriangles(wr, Rot(a), coords);
//       Wr.PutText(wr, "end \n"); 
//       Wr.PutText(wr, "object surfac scale 2 2 2 \n");       
//       Wr.Close(wr)
//     END
//   END WriteRayTrace;
//   
// PROCEDURE WriteFileEdges(t: Topology; filename: TEXT) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     WITH 
//       name = filename & ".edge",
//       wr = FileStream.OpenWrite(name)
//     DO
//       FOR i := 0 TO t.NE-1 DO 
//         Oct.PrintEdge(wr, t.edge[i].edge);
//       END;
//       Wr.Close(wr)
//     END
//   END WriteFileEdges;
//   
// PROCEDURE ReadFileEdges(t: Topology; filename: TEXT) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     WITH 
//       name = filename & ".edge",
//       rd = FileStream.OpenWrite(name)
//     DO
//       FOR i := 0 TO t.NE-1 DO 
//         Oct.ReadEdge(rd, t.edge[i].edge, map);
//       END;
//       Wr.Close(wr)
//     END
//   END ReadFileEdges;
//  
// PROCEDURE CalcTriang(a: Arc; shape: TEXT) =
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     Enum(a, GluePatch);
//     WITH 
//       s = Midle(a),
//       t = MakeTopology(s),
//       coords = NEW(REF Coords, t.NV)
//     DO
//       Enum(, prt);
//       WriteFileEdges(t, "ArquivoTeste");
//       Wr.PutText(stderr, "\n");            
// (*     SimplexOpt(s, coords^, shape, 1000); *)
// (*     GradOpt(s, coords^, shape, 1000);   *)  
//      PraxOpt(s, t, coords^, shape, 10000);   
// (*     GradTOpt(s, coords^, shape, 1000); *)         
// (*      WriteWireFrame(s, coords^, shape, "F");*)
//     END;
//   END CalcTriang;
// 
// PROCEDURE MakeTetra(): Arc =
//   VAR s, t, a: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     s := MakeRing(4,1);
//     a := MakeSOct(1);
//     t := Rnext(Rnext(s));
//     
//     Splice(s, a); Splice(t, Sym(a));
//     a := MakeSOct(1); 
//     Splice(Sym(s), a); Splice(Sym(t), Sym(a));
// 
//     CalcTriang(a, "tetra");
//     Wr.PutText(stderr, "Exit MakeTetra:\n"); 
//     RETURN a    
//   END MakeTetra;
// 
// PROCEDURE MakeCubeRamid(): Arc =
//   VAR b, t, e: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     b := MakeRing(4,1);
//     t := MakeRing(4,1); t := Oprev(t);
//     FOR i := 1 TO 4 DO
//       e := MakeSOct(1);
//       Splice(b, e); Splice(Sym(e), t);
//       b := Rprev(b); t := Rnext(t)
//     END;
//     CalcTriang(b, "ramid");
//     Wr.PutText(stderr, "Exit MakeCubeRamid:\n"); 
//     RETURN b    
//   END MakeCubeRamid;
// 
// PROCEDURE MakeStick ():Arc =
//   VAR a: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MakeSOct(GridOrder);
//     Splice(a, Sym(a));
//     CalcTriang(a, "stick");
//     Wr.PutText(Stdio.stderr, "Exit MakeStick:\n");
//     RETURN a;
//   END MakeStick;
// 
// PROCEDURE MakeRing (n: INTEGER; order: CARDINAL): Arc =
//   VAR fst, a, b: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MakeSOct(order); fst := a;
//     FOR i := 2 TO n DO
//       b := MakeSOct(order); Splice(b, Sym(a));
//       a := b;
//     END;
//     Splice(fst, Sym(a));
//     RETURN fst
//   END MakeRing;
// 
// PROCEDURE MakeCube (): Arc =
//   VAR b, t, e: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     b := MakeRing(4,GridOrder);
//     t := MakeRing(4,GridOrder); t := Oprev(t);
//     FOR i := 1 TO 4 DO
//       e := MakeSOct(GridOrder);
//       Splice(b, e); Splice(Sym(e), t);
//       b := Rprev(b); t := Rnext(t)
//     END;
//     CalcTriang(b, "cube");
//     Wr.PutText(stderr, "Exit MakeCube:\n"); 
//     RETURN b    
//   END MakeCube;
// 
// PROCEDURE MakeSausage (): Arc =
//   VAR s, t, e: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     t := MakeRing(2,GridOrder);
//     FOR i := 1 TO 2 DO
//       t := Oprev(t);
//       s := MakeRing(2,GridOrder);
//       FOR j  := 1 TO 2 DO 
//         e := MakeSOct(GridOrder);
//         Splice(t, e); Splice(Sym(e), s);
//         t := Lnext(t); s := Rnext(s)
//       END;
//       t := s;
//     END;
//     CalcTriang(s, "sausage");
//     Wr.PutText(stderr, "Exit sausage:\n"); 
//     RETURN s    
//   END MakeSausage;
// 
// PROCEDURE MakePPlane (): Arc =
//   VAR e: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *> 
//   BEGIN
//     e := MakeSOct(GridOrder);
//     Splice (Flip(Sym(e)), e);
//     CalcTriang(e, "pplane");
//     Wr.PutText(stderr, "Exit MakePPlane: \n");
//     RETURN e
//   END MakePPlane;
// 
// PROCEDURE MakeTorus(): Arc =
//   VAR a, b: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MakeSOct(GridOrder);
//     b := MakeSOct(GridOrder);
//     Splice(a, b);
//     Splice(Sym(a), a);
//     Splice(Sym(b), a);
//     CalcTriang(a, "torus");
//     Wr.PutText(stderr, "Exit MakeTorus: \n");    
//     RETURN a
//   END MakeTorus;
//  
// PROCEDURE MakeKlein(): Arc =
//   VAR a, b: Arc;
//   <* FATAL Wr.Failure, Thread.Alerted *>
//   BEGIN
//     a := MakeSOct(GridOrder);
//     b := MakeSOct(GridOrder);
//     Splice(a, b);
//     Splice(Sym(a), a);
//     Splice(Flip(Sym(b)), a);
//     CalcTriang(a, "klein");
//     Wr.PutText(stderr, "Exit MakeKlein: \n");
//     RETURN a
//   END MakeKlein;
//   
// PROCEDURE Main() =  
//   <* FATAL Wr.Failure, Thread.Alerted *>  
//   BEGIN
//   (*  
//      cube := MakeCube();
//      FOR i := 1 TO 4 DO
//        prt(cube); cube := Flip(cube);
//        prt(cube); cube := Flip(cube);
//        cube := Rot(cube)
//       END; 
//   *)
//      tetra := MakeTetra();   
//     (* ramid := MakeCubeRamid(); *)  
//     (* sausage := MakeSausage();  *)  
//     (* stick := MakeStick();  *)
//     (* pplane := MakePPlane(); *)
//     (* ring := MakeRing(4,GridOrder);  *)
//     (* klein := MakeKlein();  *)
//     (* torus := MakeTorus(); *)   
//   END Main;
// 
// BEGIN 
//   Main();
// END TestOctet.

#endif

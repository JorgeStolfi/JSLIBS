(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.ansp.br>    *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.ansp.br>  *)
(*   Jorge Stolfi        - DEC Systems Research Center <stolfi@src.dec.com> *)
(*                                                                          *)
(* This file can be freely distributed, modified, and used for any          *)
(*   non-commercial purpose, provided that this copyright and authorship    *)
(*   notice be included in any copy or derived version of this file.        *)
(*                                                                          *)
(* DISCLAIMER: This software is offered ``as is'', without any guarantee    *)
(*   as to fitness for any particular purpose.  Neither the copyright       *)
(*   holder nor the authors or their employers can be held responsible for  *)
(*   any damages that may result from its use.                              *)
(****************************************************************************)

(* Last modified on Thu Aug 20 16:39:39 PDT 1992 by stolfi                  *)

INTERFACE ShrinkDAG;

(*
  Exports procedures to analyze the pointer distribution of a DAG.T,
  and to compress/uncompress it using various pointer encoding tricks.
  *)

IMPORT DAG, Wr, Encoding;
FROM Basics IMPORT NAT, Letter;

TYPE 
  BYTE = BITS 8 FOR [0..255];
  
  PtrField = {Dest, Rest};      (* Pointer fields of a DAG node *)
  PtrCounts = ARRAY [0..256*256] OF NAT;

  PtrCoding = RECORD
      max_val: NAT;    (* Maximum pointer value *)
      pop_val: NAT;    (* Most popular pointer value *)
      pop_mask: BYTE;  (* Bit of /rd/ byte to use for /pop_val/, or 0 if not avail. *)
      n1: NAT;         (* Pointer values [0..n1-1] should be coded with 1 byte *)
      n2: NAT;         (* Pointer values [n1..n1+n2-1] should be coded with 2 bytes *)
      tot_bytes: NAT;  (* Total size (bytes) of pointers after encoding. *)
    END;
    
  LetField = {Rd, Wr};                 (* Letter fields of a DAG node *)
  LetMap = ARRAY Letter OF Letter;     (* Compact coding of each Letter *)
  LetCounts = ARRAY Letter OF NAT;
  
  LetCoding = RECORD (* Frequencies for the /rd/ or /wr/ letter fields *)
      num: NAT;             (* Number of distinct letters *)
      map: REF LetMap;      (* Compact recoding of each letter *)
    END;
    
  DAGCoding = RECORD
      nodes: NAT;                        (* Total number of nodes *)
      let: ARRAY LetField OF LetCoding;  (* Statistics for /rd/ and /wr/ fields *)
      ptr: ARRAY PtrField OF PtrCoding;  (* Statistics for /dest/ and /rest/ *)
      tot_bytes: NAT;                    (* Compressed DAG size (excluding header) *)
    END;

  DAGStats = RECORD
      let_ct: ARRAY LetField OF REF LetCounts; (* Frequency of each letter *)
      ptr_ct: ARRAY PtrField OF REF PtrCounts; (* Frequency of each pointer value *)
      coding: DAGCoding;                       (* Best coding given the observed stats. *)
    END;

PROCEDURE Analyze(dag: DAG.T): DAGStats;
(*
  Returns a table /ds.ptr/ with statistics for the /rd/ and /wr/
  letters and the /rest/ and /dest/ pointers in /dag/.
  
  Specifically, /ds.ptr[f].ct[v]/ tells how many /f/ pointers
  have "value" /v/, where the value is defined by the coding
  function.  The entry /ct[v]/ for /v=65536/ includes also all pointers with
  value greater than /65536/.
  
  In addition to the raw counts, the result /ds/ also contains 
  a set of parameters /ds.coding/ needed by the encoding and decoding procedures.
  For each letter field /f/, the record /ds.coding.let[f]/
  tells the number /num/ of distinct letter values found, and
  a table /map/ that maps them to the range [0..num-1].
  For each pointer field /f/, the record /ds.coding.ptr[f]/ contains:
|  
|     /max_val/ = maximum observed pointer value
|
|     /pop_val/ = most popular pointer value in [0..65535]
|
|     /rd_mark/ = mask for the bit of the /rd/ letter that is to be used
|                 to encode /pop_val/; or 0 if such bit is not available.
|
|     /n1, n2/ = number of pointer values to be encoded in 1 and 2 bytes, resp.
|
|     /tot_bytes/ = total bytes needed to encode these pointer values,
|             for the optimal choice of /s/ and /m/.
|
  *)
  
PROCEDURE PrintStats(wr: Wr.T; ds: DAGStats; rd_e, wr_e: Encoding.T := NIL);
(*
  Prints a summary of the statistics in /ds/ to the given writer.
  The encodings /rd_e/ and /wr_e/ are used only to print the letters legibly, 
  and default to PlainEncoding.New(). *)
  
PROCEDURE Encode(dag: DAG.T; coding: DAGCoding): REF ARRAY OF BYTE;
(*
  Returns a compact binary representation of /dag/. *)
  
PROCEDURE Decode(READONLY b: ARRAY OF BYTE; coding: DAGCoding): DAG.T;
(*
  Expands the given byte encoding of a DAG.T (as created by Compress)
  and returns the corresponding DAG.T. *)

END ShrinkDAG.




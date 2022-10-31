#ifndef dnae_seq_H
#define dnae_seq_H

/* Numerically encoded and filtered DNA sequences */
/* Last edited on 2022-10-31 09:38:20 by stolfi */

#define dnae_seq_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <vec.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>

#include <dnae_sample.h>
#include <dnae_datum.h>

/* This interface defines a discrete representation for 
  signals, as strings of /datum/ elements ({dnae_datum_t}s). */

typedef msm_seq_id_t dnae_seq_id_t; 
  /* A {dnae_seq_id_t} is an internal identifier for a sequence, 
     presently an index into a table of sequences. */

#define dnae_seq_id_none msm_seq_id_none
  /* A {dnae_seq_id_t} value that means "no sequence". */

typedef struct dnae_seq_t
  { msm_seq_desc_t sd;        /* Abstract sequence attributes for matching. */
    char *cmt;                /* Comment text. */
    dnae_datum_scale_t sfac;  /* Per-channel scale factors for decoded samples. */
    dnae_datum_vec_t dv;      /* The datum array. */
  } dnae_seq_t;  
  /* A multichannel signal, sampled at and encoded. */

/* CREATION */

dnae_seq_t dnae_seq_new(int32_t n);
  /* Allocates a new datum sequence with {n} sample datums, initially
    all set to {(0,..0)}. The sequence will have {sd.id=dnae_seq_id_none},
    {sd.name=cmt=NULL}, {sd.rev=NULL}, {sd.estep=0}, {sd.skip=0}, {sd.size=n}, and
    {sfac[0..dnae_CHANNELS-1]=1}. */

dnae_seq_t dnae_seq_from_datum_vec
  ( dnae_seq_id_t id, 
    char *name, 
    bool_t rev,
    int8_t estep, 
    int32_t skip, 
    char *cmt, 
    dnae_datum_scale_t *sfac,
    dnae_datum_vec_t dv
  );
  /* Assembles a sequence descriptor {sd} with the given fields.
    The number of sequence positions {sd.size} is the number of datums in {dv}. */

dnae_seq_t dnae_seq_from_nucleic_string(dnae_seq_id_t id, char *name, bool_t rev, char *cmt, char *bas);
  /* Converts the DNA/RNA string {bas} into a {dnae_datum_t} sequence.
    Each base will become one datum, with no filtering or resampling.
    If {rev} is true, applies reversal-complement to it first.
    The sequence's numeric and alphabetic identifiers will be {id} and
    {name}, respectively. The comment text will be {txt}. */

dnae_seq_t dnae_seq_copy(dnae_seq_t *seq);
  /* Creates a heap copy of {*seq}, including all its internal storage
    (so that {dnae_seq_free_datums(seq)} will not affect the copy, and vice-versa). */

dnae_seq_t dnae_seq_copy_sub(dnae_seq_t* seq, int32_t ini, int32_t fin);
 /* Creates a heap copy of the segment of {*seq} from datum number {ini} to 
    datum {fin} inclusive. The datum vector is newly allocated 
    so that {dnae_seq_free_datums(seq)} will not affect the copy, and vice-versa. */

/* EXTRACTING DATA */

int32_t dnae_seq_num_datums(dnae_seq_t *seq);
  /* Number of datums in sequence {seq}. */

dnae_sample_enc_t *dnae_seq_get_sample_enc_address(dnae_seq_t *seq, int32_t i, int32_t k);
dnae_sample_enc_t dnae_seq_get_sample_enc(dnae_seq_t *seq, int32_t i, int32_t k);
void dnae_seq_set_sample_enc(dnae_seq_t *seq, int32_t i, int32_t k, dnae_sample_enc_t s);
  /* These procedures return the address, return the encoded value, and set the
    encoded value of the sample in channel {k} of datum {i} from
    the sequence {*seq}.  The index {i} must lie in the range
    {[0..seq->dv.ne-1]}. */
 
dnae_datum_t *dnae_seq_get_datum_address(dnae_seq_t *seq, int32_t i);
dnae_datum_t dnae_seq_get_datum(dnae_seq_t *seq, int32_t i);
void dnae_seq_set_datum(dnae_seq_t *seq, int32_t i, dnae_datum_t d);
  /* These procedures return the address, return the value, and set
    the value of encoded datum {i} from the sequence {*seq}. The index {i}
    must lie in the range {[0..seq->dv.ne-1]}. */

dnae_datum_t dnae_seq_eval(dnae_seq_t *seq, double x, bool_t smooth);
  /* Evaluates the sequence {seq} for argument {x}, which must be
    in the range {[0 _ seq.size-1]}.  If {x} is not an
    integer, performs linear interpolation of
    neighboring decoded datums, and encodes the result using the
    same scale factors {seq->sfac}. */

/* SEQUENCE I/O

  The procedures in this section read or write a datum sequence {seq} from
  a file. The file format is described by the string {dnae_seq_file_format_INFO} below.
  */
  
#define dnae_seq_type_name "encoded_bio_seq"
#define dnae_seq_version "2014-06-14"
  /* Current format, without {circular}, with {reverse}, 2-byte encoded samples. */

#define dnae_seq_version_old_2 "2013-10-23"
  /* Old format, without {circular}, with {reverse}, 1-byte encoded samples. */

#define dnae_seq_version_old_1 "2008-01-29" 
  /* Old format, with {circular} and no {reverse}, 1-byte encoded samples. */

#define dnae_seq_file_format_INFO \
  "A sequence file consists of \n" \
  "\n" \
  "      a standard header line, \"begin {tname} (version of {vdate})\", where" \
  " {tname} is {dnae_seq_type_name = \"" dnae_seq_type_name "\"}" \
  " and {vdate} is {dnae_seq_version = \"" dnae_seq_version "\"};" \
  "\n" \
  "      zero or more comment lines, starting with \"|\";\n" \
  "\n" \
  "      a set of sequence parameter lines:\n" \
  "        \"id = {NUMID}\"\n"\
  "        \"name = {NAME}\"\n" \
  "        \"reversed = {REVERSED}\"\n" \
  "        \"resampling = {RESAMPLING}\"\n" \
  "        \"offset = {OFFSET}\"\n" \
  "        \"channels = {NC}\"\n" \
  "        \"scale = {SFAC[0]} {SFAC[1]} .. SFAC[NC-1]}\"\n" \
  "        \"datums = {NS}\"\n" \
  "\n" \
  "      then {NS} datum lines, each with {NC} samples"\
  " in the format \"{SMP[i,0]} {SMP[i,1] .. {SMP[i,NC-1]}\"\n" \
  "      for {i} ranging from 0 to {NS-1};\n" \
  "\n" \
  "      a standard footer line, \"end {tname}\".\n" \
  "\n" \
  "  In the paramter lines," \
  "\n" \
  "      {NUMID} (field {seq.sd.id}) is an ID number"\
  " of the original biosequence from which this datum"\
  " sequence is derived, e.g. an index into some"\
  " sequence database.\n"\
  "\n" \
  "      {NAME} (field {seq.sd.name}) is a menemonic tag"\
  " for that sequence, without quotes or embedded blanks.\n"\
  "\n" \
  "      {REVERSED} (field {seq.sd.rev} is a logical"\
  " value \"F\" or \"T\" (without quotes), specifying"\
  " whether this datum sequence is derived from the original"\
  " biosequence directly, or after reversing and"\
  " complementing it.\n"\
   "\n" \
 "       {RESAMPLING} (field {seq.sd.estep}) is a signed integer"\
  " that defines the total amount of subsampling/supersampling"\
  " that was applied to the original biosequence.  Specifically,"\
  " each step of this sequence corresponds to {2^RESAMPLING}"\
  " steps of the original.  The value is negative for"\
  " supersampling, positive for downsampling.\n"\
  "\n" \
  "      {OFFSET} (field {seq.sd.skip}) is a signed integer that"\
  " specifies the position of the first datum relative to the first"\
  " base of the original biosequence.  Specifically, datum {i} of"\
  " this sequence corresponds to element {(i - skip)*2^RESAMPLING} of"\
  " the original biosequence.  Note that if {RESAMPLING} is negative"\
  " the corresponding index may be fractionary.\n"\
  "\n" \
  "      {NC} is a positive integer, the number of samples per datum.\n"\
  "\n" \
  "      {SFAC[k]} is a real scale factor to be multiplied into all"\
  " samples of the corresponding" \
  " channel, after applying {dnae_sample_decode}.\n"\
  "\n" \
  "      {NS} (field {seq.sd.size) is a positive integer, the number of datums"\
  " in this sequence.\n"\
  "\n" \
  "      {SMP[i,k]} is a small signed integer, the un-scaled and encoded"\
  " sample value in channel {k} of datum {i}.\n" \
  "\n" \
  "  The format changed in "dnae_seq_version "; the previous"\
  " versions (" dnae_seq_version_old_1 ", " dnae_seq_version_old_1 ") had"\
  " different header format and 1-byte (instead of 2-byte) encoded samples.\n" \
  "\n" \
  ""     

dnae_seq_t dnae_seq_read(FILE *rd);
  /* Reads a numeric DNA/RNA sequence from file {rd},
    and turns it into a sequence. The sequence will have 
    {id=dnae_seq_id_none}, {name=NULL}. */
  
dnae_seq_t dnae_seq_read_named(char *fname, char *tag, char *ext);
  /* Same as {dnae_seq_read}, but reads from a file named "{fname}{tag}{ext}".  If that
    file name is "-", however, reads from {stdin}. */
  
void dnae_seq_write(FILE *wr, dnae_seq_t *seq);
  /* Writes a numeric DNA/RNA sequence to file {wr}. */

void dnae_seq_write_named(char *fname, char *tag, char *ext, dnae_seq_t *seq);
  /* Same as {dnae_seq_write}, but creates a disk file called "{name}{tag}{ext}".
    If that file name is "-", however, writes to {stdout}. */

/* DNA/RNA SEQUENCE I/O */

dnae_seq_t dnae_seq_read_from_nucleic_file(FILE *rd, dnae_seq_id_t id, char *name, bool_t rev);
 /* Reads a DNA/RNA sequence from {rd} and converts it to a numerically
    encoded {dnae_seq_t}, with no filtering or resampling. The sequence's numeric and 
    alphabetic identifiers will be {id} and {name}, respectively. If {rev} is true, 
    the sequence is reversed and complemented before encoding. */

dnae_seq_t dnae_seq_read_from_nucleic_file_named(char *fname, char *tag, char *ext, dnae_seq_id_t id, char *name, bool_t rev);
  /* Same as {dnae_seq_read_from_nucleic_file}, but reads from disk file "{fname}{tag}{ext}"  */
  
void dnae_seq_free_datums(dnae_seq_t *seq);
   /* Reclaims the internal storage of sequence {*seq} (but not {*seq} itself). */

void dnae_seq_free(dnae_seq_t *seq);
   /* Reclaims the sequence record {*seq} and all its internal storage. */
     
void dnae_seq_multi_free_datums(dnae_seq_t seq[], int32_t maxLevel);
  /* Reclaims the internal storage of sequences {seq[0..maxLevel]} 
    (but not the array {*seq} itself). */
   
void dnae_seq_multi_free(dnae_seq_t *seq[], int32_t maxLevel);
  /* Reclaims the sequences {*(seq[0..maxLevel])} and all
    their internal storage. */
    
/* FILTERING AND INTERPOLATING */

dnae_seq_t dnae_seq_filter(dnae_seq_t *seq, double_vec_t *wtb, int8_t ek, char *wcmt);
  /* Returns a sequence {new} obtained by filtering a numeric datum sequence {seq} 
    with the kernel {wtb}, and resampling it with step size {2^ek} times its
    current step.  The exponent {ek} must be non-negative.
    
    The descriptor of {new} is given by {msm_desc_filter(seq.sd,nw,ek)}
    where {nw=wtb->ne}.
    
    More precisely the datums {pnew[0..nnew-1]} of the new sequence
    {new} are related to the datums {pseq[0..nseq-1]} of the old ones by
    the formula
    
      {pnew[i] = SUM{ pseq[2^ek*i+j+d] * wtb[j] : j in 0..nw-1}}
      
    where {d=new.skip*2^ek-hw} and {hw = (nw-1)/2}.
    
    The string {wcmt}, if not NULL, is assumed to be the name or
    description of the weight table {wtb}, and is appended to the
    comment {seq->cmt} to make the new string's comment. For each
    channel {k}, the new scale factor {sfac.f[k]} will be the
    root-mean-square value of the samples in channel {k} of the filtered
    sequence. */
    
dnae_seq_t dnae_seq_interpolate(dnae_seq_t *seq, int8_t ek);
  /* Returns a sequence {new} obtained by interpolating a numeric datum sequence {seq} 
    with step size {2^ek} times its current step.  The exponent {ek} must be zero or negative.
    
    The descriptor of {new} is given by {msm_seq_desc_resample(seq.sd,ek)}.
    
    More precisely, if {pnew[0..nnew-1]} are the datums of the
    new sequence {new}, and {pseq[o..nseq-1]} are the datums of {seq},
    then {pnew[i*stepn] == pseq[i]} for {i} in {0..nseq-1}.
    The other elements of {pnew} are interpolated from the 
    adjacent ones by an interpolating C1 cubic spline.
    
    Also appends the text "interpolated with step {1/stpn}" to the
    comment {seq->cmt} to make the new string's comment. For each
    channel {k}, the new scale factor {sfac.f[k]} will be the
    root-mean-square value of the samples in channel {k} of the
    interpolated sequence.
    
    Currently uses C1 cubic Hermite interpolation in each interval,
    with slopes at the original elements estimated by finite differences. */
    
#endif

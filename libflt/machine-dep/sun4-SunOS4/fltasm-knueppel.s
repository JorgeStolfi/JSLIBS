! Setting the IEEE rounding mode flags on a Sun SPARCstation
! See flt.h
! Created by J. Stolfi on 93-12-21
! Based on fpSparc.s version 1.2 of 93-14-06 by O. Knueppel
! Careful: these routines do not update _fp_direction, which
!   is sinful according to Douglas Priest of Sun.

! The following should go in flt.h:
! /* Setting the IEEE rounding mode bits on a SPARC: */
! 
! #define ROUND_DOWN  flt_round_down()
! #define ROUND_UP    flt_round_up()
! #define ROUND_NEAR  flt_round_near()
! #define ROUND_ZERO  flt_round_zero()
! 
! void flt_round_down (void);
! void flt_round_up   (void);
! void flt_round_near (void);
! void flt_round_zero (void);
!   /* These routines set the IEEE FP rounding direction. */
!   /* They should be a lot faster than "ieee_flags". */

      .data
      .align  4
_flt_sr_round_near:
      .word   0x00000000
_flt_sr_round_zero:
      .word   0x40000000
_flt_sr_round_up:
      .word   0x80000000
_flt_sr_round_down:
      .word   0xC0000000
      .text

      .align  4
      .global _flt_round_near
      .proc   020
_flt_round_near:
      sethi   %hi(_flt_sr_round_near),%g2
      retl
      ld      [%g2+%lo(_flt_sr_round_near)],%fsr

      .align  4
      .global _flt_round_zero
      .proc   020
_flt_round_zero:
      sethi   %hi(_flt_sr_round_near),%g2
      retl
      ld      [%g2+%lo(_flt_sr_round_near)],%fsr

      .align  4
      .global _flt_round_up
      .proc   020
_flt_round_up:
      sethi   %hi(_flt_sr_round_up),%g2
      retl
      ld      [%g2+%lo(_flt_sr_round_up)],%fsr

      .align  4
      .global _flt_round_down
      .proc   020
_flt_round_down:
      sethi   %hi(_flt_sr_round_down),%g2
      retl
      ld      [%g2+%lo(_flt_sr_round_down)],%fsr


      .align 4
      .global _flt_get_fsr
      .proc 012
_flt_get_fsr: 
      !#PROLOGUE# 0
      save %sp,-120,%sp
      !#PROLOGUE# 1
      st %fsr,[%fp-20]
      ld [%fp-20],%i0
      ret
      restore

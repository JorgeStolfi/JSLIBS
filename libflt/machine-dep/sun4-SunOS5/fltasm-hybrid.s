! Setting the IEEE rounding mode bits on a SPARC
! Jorge Stolfi 93-12-26

! The following declarations should go in flt.h:	

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
! 
! int flt_get_fsr(void);
!   /* Returns the FP status register. (For debugging) */

      .data
      .align 4
_flt_sr_round_near:
      .word 0x00000000
_flt_sr_round_zero:
      .word 0x40000000
_flt_sr_round_up:
      .word 0x80000000
_flt_sr_round_down:
      .word 0xC0000000
      .text

      .align 4
      .global flt_round_near
      .proc 020
flt_round_near:
      set 0,%o0
      sethi %hi(fp_direction),%g2
      st %o0,[%g2+%lo(fp_direction)]
      sethi %hi(_flt_sr_round_near),%g2
      retl
      ld [%g2+%lo(_flt_sr_round_near)],%fsr

      .align 4
      .global flt_round_zero
      .proc 020
flt_round_zero:
      set 1,%o0
      sethi %hi(fp_direction),%g2
      st %o0,[%g2+%lo(fp_direction)]
      sethi %hi(_flt_sr_round_near),%g2
      retl
      ld [%g2+%lo(_flt_sr_round_near)],%fsr

      .align 4
      .global flt_round_up
      .proc 020
flt_round_up:
      set 2,%o0
      sethi %hi(fp_direction),%g2
      st %o0,[%g2+%lo(fp_direction)]
      sethi %hi(_flt_sr_round_up),%g2
      retl
      ld [%g2+%lo(_flt_sr_round_up)],%fsr

      .align 4
      .global flt_round_down
      .proc 020
flt_round_down:
      set 3,%o0
      sethi %hi(fp_direction),%g2
      st %o0,[%g2+%lo(fp_direction)]
      sethi %hi(_flt_sr_round_down),%g2
      retl
      ld [%g2+%lo(_flt_sr_round_down)],%fsr


      .align 4
      .global flt_get_fsr
      .proc 012
flt_get_fsr: 
      !#PROLOGUE# 0
      save %sp,-120,%sp
      !#PROLOGUE# 1
      st %fsr,[%fp-20]
      ld [%fp-20],%i0
      ret
      restore

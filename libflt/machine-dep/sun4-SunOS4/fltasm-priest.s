! Setting the IEEE rounding mode bits on a SPARC
! Jorge Stolfi 93-12-26
! Based on suggestions by Douglas Priest

! The following should go in flt.h:
! /* Setting the IEEE rounding mode bits on a SPARC: */
! 
! #define ROUND_DOWN  flt_set_rounding (fp_negative)
! #define ROUND_UP    flt_set_rounding (fp_positive)
! #define ROUND_NEAR  flt_set_rounding (fp_nearest)
! #define ROUND_ZERO  flt_set_rounding (fp_tozero)
! 
! void flt_set_rounding (enum fp_direction_type dir);
!   /* This routine sets the IEEE FP rounding direction. */
!   /* It should be faster than "ieee_flags". */

	.text
	.align 4
	.global _flt_set_rounding
_flt_set_rounding:
	!#PROLOGUE# 0
	!#PROLOGUE# 1
	and %o0,3,%o0      ! Ensure argument is in [0..3]
	sethi %hi(_fp_direction),%o1     ! Update fp_direction
	st %o0,[%o1+%lo(_fp_direction)]  ! "
	sll %o0,30,%o0     ! Shift new rounding dir to proper place
	set 0xc0000000,%o1 ! Mask of rounding direction bits
        st %fsr,[%sp+0x44]   ! Get FP status register
	ld [%sp+0x44],%o2    !  (must we go through memory?)
	andn %o2,%o1,%o2
        or %o0,%o2,%o0
	st %o0,[%sp+0x44]    ! Set FP status register
	ld [%sp+0x44],%fsr   !  (must we go through memory?)
	retl
	nop


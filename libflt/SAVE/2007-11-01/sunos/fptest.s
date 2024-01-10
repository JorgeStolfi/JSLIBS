.text
	.align 4
	.global set_rounding
	.proc	020
set_rounding:
	!#PROLOGUE# 0
	save %sp,-120,%sp
	!#PROLOGUE# 1
	and %i0,3,%o0      ! Ensure argument is in [0..3]
	sll %o0,30,%o0     ! Shift new rounding dir to proper place
	set 0xc0000000,%o1 ! Mask of rounding direction bits
        st %fsr,[%fp+68]   ! Get FP status register
	ld [%fp+68],%o2    !  (must we go through memory?)
	andn %o2,%o1,%o2
        or %o0,%o2,%o0
	st %o0,[%fp+68]    ! Set FP status register
	ld [%fp+68],%fsr   !  (must we go through memory?)
	ret
	restore

	.align 4
	.global get_fsr
	.proc	012
get_fsr:	
	!#PROLOGUE# 0
	save %sp,-120,%sp
	!#PROLOGUE# 1
	st %fsr,[%fp-20]
	ld [%fp-20],%i0
	ret
	restore

	

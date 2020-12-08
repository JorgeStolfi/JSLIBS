dspmat_extra.ho: dspmat_extra.h
dspmat_extra.ho: /home/jstolfi/include/bool.ho
dspmat_extra.ho: dspmat.ho

spmat_io_def.ho: spmat_io_def.h
spmat_io_def.ho: /home/jstolfi/include/affirm.ho
spmat_io_def.ho: /home/jstolfi/include/fget.ho
spmat_io_def.ho: /home/jstolfi/include/nget.ho
spmat_io_def.ho: /home/jstolfi/include/filefmt.ho
spmat_io_def.ho: spmat.ho
spmat_io_def.ho: spmat_io.ho

dspmat.ho: dspmat.h
dspmat.ho: /home/jstolfi/include/bool.ho
dspmat.ho: spmat.ho
dspmat.ho: spmat_io.ho
dspmat.ho: spmat_linalg.ho

dspmat_linsys_GS.ho: dspmat_linsys_GS.h
dspmat_linsys_GS.ho: /home/jstolfi/include/bool.ho
dspmat_linsys_GS.ho: dspmat.ho

spmat_linalg_def.ho: spmat_linalg_def.h
spmat_linalg_def.ho: spmat.ho
spmat_linalg_def.ho: spmat_linalg.ho

spmat.ho: spmat.h
spmat.ho: spmat.ho
spmat.ho: spmat.ho
spmat.ho: spmat_def.ho

spmat_io.ho: spmat_io.h
spmat_io.ho: spmat_io.ho
spmat_io.ho: spmat_io.ho
spmat_io.ho: spmat.ho
spmat_io.ho: spmat_io_def.ho

spmat_linalg.ho: spmat_linalg.h
spmat_linalg.ho: spmat_linalg.ho
spmat_linalg.ho: spmat_linalg.ho
spmat_linalg.ho: spmat.ho
spmat_linalg.ho: spmat_linalg_def.ho

dspmat_linsys_ALT.ho: dspmat_linsys_ALT.h
dspmat_linsys_ALT.ho: /home/jstolfi/include/bool.ho
dspmat_linsys_ALT.ho: dspmat.ho

spmat_def.ho: spmat_def.h
spmat_def.ho: spmat.ho

dspmat_extra.o: dspmat_extra.c
dspmat_extra.o: /home/jstolfi/include/bool.ho
dspmat_extra.o: /home/jstolfi/include/jsmath.ho
dspmat_extra.o: dspmat.ho
dspmat_extra.o: dspmat_linsys_GS.ho
dspmat_extra.o: dspmat_linsys_ALT.ho
dspmat_extra.o: dspmat_extra.ho

spmat_io.o: spmat_io.c
spmat_io.o: /home/jstolfi/include/affirm.ho
spmat_io.o: /home/jstolfi/include/filefmt.ho
spmat_io.o: /home/jstolfi/include/nget.ho
spmat_io.o: /home/jstolfi/include/fget.ho
spmat_io.o: spmat.ho
spmat_io.o: spmat_io.ho

dspmat.o: dspmat.c
dspmat.o: /home/jstolfi/include/bool.ho
dspmat.o: spmat.ho
dspmat.o: spmat_io.ho
dspmat.o: spmat_linalg.ho
dspmat.o: dspmat.ho

dspmat_linsys_GS.o: dspmat_linsys_GS.c
dspmat_linsys_GS.o: /home/jstolfi/include/bool.ho
dspmat_linsys_GS.o: /home/jstolfi/include/affirm.ho
dspmat_linsys_GS.o: /home/jstolfi/include/jsmath.ho
dspmat_linsys_GS.o: dspmat.ho
dspmat_linsys_GS.o: dspmat_linsys_GS.ho

spmat.o: spmat.c
spmat.o: spmat.ho
spmat.o: /home/jstolfi/include/affirm.ho

spmat_linalg.o: spmat_linalg.c
spmat_linalg.o: spmat.ho
spmat_linalg.o: spmat_linalg.ho
spmat_linalg.o: /home/jstolfi/include/affirm.ho

dspmat_linsys_ALT.o: dspmat_linsys_ALT.c
dspmat_linsys_ALT.o: /home/jstolfi/include/bool.ho
dspmat_linsys_ALT.o: /home/jstolfi/include/affirm.ho
dspmat_linsys_ALT.o: dspmat.ho
dspmat_linsys_ALT.o: /home/jstolfi/include/jsmath.ho
dspmat_linsys_ALT.o: dspmat_linsys_ALT.ho


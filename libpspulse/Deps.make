psp_pulse_mother.ho: psp_pulse_mother.h
psp_pulse_mother.ho: psp_basic.ho
psp_pulse_mother.ho: psp_grid.ho
psp_pulse_mother.ho: psp_pulse.ho

psp_pulse_mother_B.ho: psp_pulse_mother_B.h
psp_pulse_mother_B.ho: psp_pulse.ho

psp_basic.ho: psp_basic.h
psp_basic.ho: spline_basic.ho
psp_basic.ho: rgrid_basic.ho
psp_basic.ho: /home/jstolfi/include/bz_basic.ho

psp_tent.ho: psp_tent.h
psp_tent.ho: psp_pulse.ho
psp_tent.ho: psp_basic.ho
psp_tent.ho: /home/jstolfi/include/vec.ho
psp_tent.ho: /home/jstolfi/include/bool.ho

psp_pulse.ho: psp_pulse.h
psp_pulse.ho: psp_basic.ho
psp_pulse.ho: psp_grid.ho

psp_tent_mother.ho: psp_tent_mother.h
psp_tent_mother.ho: psp_pulse.ho
psp_tent_mother.ho: psp_basic.ho
psp_tent_mother.ho: /home/jstolfi/include/vec.ho
psp_tent_mother.ho: /home/jstolfi/include/bool.ho

psp_pulse_mother_H.ho: psp_pulse_mother_H.h
psp_pulse_mother_H.ho: psp_pulse.ho

spline_basic.ho: spline_basic.h

rgrid_basic.ho: rgrid_basic.h
rgrid_basic.ho: spline_basic.ho
rgrid_basic.ho: /home/jstolfi/include/interval.ho
rgrid_basic.ho: /home/jstolfi/include/box.ho
rgrid_basic.ho: /home/jstolfi/include/bool.ho
rgrid_basic.ho: /home/jstolfi/include/sign.ho
rgrid_basic.ho: /home/jstolfi/include/vec.ho
rgrid_basic.ho: /home/jstolfi/include/indexing.ho

psp_pulse_mother_N.ho: psp_pulse_mother_N.h
psp_pulse_mother_N.ho: psp_basic.ho
psp_pulse_mother_N.ho: psp_pulse.ho

psp_grid.ho: psp_grid.h
psp_grid.ho: /home/jstolfi/include/interval.ho
psp_grid.ho: /home/jstolfi/include/box.ho
psp_grid.ho: rgrid_basic.ho
psp_grid.ho: /home/jstolfi/include/bz_basic.ho
psp_grid.ho: /home/jstolfi/include/bool.ho
psp_grid.ho: /home/jstolfi/include/sign.ho
psp_grid.ho: /home/jstolfi/include/vec.ho

psp_pulse_mother.o: psp_pulse_mother.c
psp_pulse_mother.o: /home/jstolfi/include/affirm.ho
psp_pulse_mother.o: /home/jstolfi/include/bz_basic.ho
psp_pulse_mother.o: psp_basic.ho
psp_pulse_mother.o: psp_pulse_mother.ho
psp_pulse_mother.o: psp_pulse_mother_B.ho
psp_pulse_mother.o: psp_pulse_mother_H.ho
psp_pulse_mother.o: psp_pulse_mother_N.ho

psp_pulse_mother_H.o: psp_pulse_mother_H.c
psp_pulse_mother_H.o: /home/jstolfi/include/affirm.ho
psp_pulse_mother_H.o: psp_basic.ho
psp_pulse_mother_H.o: psp_pulse.ho
psp_pulse_mother_H.o: psp_pulse_mother_H.ho

psp_pulse_mother_B.o: psp_pulse_mother_B.c
psp_pulse_mother_B.o: /home/jstolfi/include/affirm.ho
psp_pulse_mother_B.o: psp_basic.ho
psp_pulse_mother_B.o: psp_pulse.ho
psp_pulse_mother_B.o: psp_pulse_mother_B.ho

psp_pulse.o: psp_pulse.c
psp_pulse.o: /home/jstolfi/include/affirm.ho
psp_pulse.o: psp_basic.ho
psp_pulse.o: psp_grid.ho
psp_pulse.o: psp_pulse_mother.ho
psp_pulse.o: psp_pulse.ho

psp_tent.o: psp_tent.c
psp_tent.o: /home/jstolfi/include/affirm.ho
psp_tent.o: /home/jstolfi/include/vec.ho
psp_tent.o: psp_basic.ho
psp_tent.o: psp_pulse.ho
psp_tent.o: psp_tent.ho

psp_pulse_mother_N.o: psp_pulse_mother_N.c
psp_pulse_mother_N.o: /home/jstolfi/include/affirm.ho
psp_pulse_mother_N.o: psp_basic.ho
psp_pulse_mother_N.o: psp_pulse.ho
psp_pulse_mother_N.o: psp_pulse_mother_N.ho


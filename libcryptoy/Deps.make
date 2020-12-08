cryptoy_file.ho: cryptoy_file.h

cryptoy_mem.ho: cryptoy_mem.h

cryptoy_file.o: cryptoy_file.c
cryptoy_file.o: cryptoy_mem.ho
cryptoy_file.o: cryptoy_file.ho

cryptoy_mem.o: cryptoy_mem.c
cryptoy_mem.o: cryptoy_mem.ho


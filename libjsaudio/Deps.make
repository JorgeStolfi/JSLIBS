jsaudio.ho: jsaudio.h

jsaudio_au.ho: jsaudio_au.h
jsaudio_au.ho: jsaudio.ho

jsaudio_dvf.ho: jsaudio_dvf.h
jsaudio_dvf.ho: jsaudio.ho

jsaudio_io.ho: jsaudio_io.h

jsaudio.o: jsaudio.c
jsaudio.o: jsaudio.ho

jsaudio_au.o: jsaudio_au.c
jsaudio_au.o: jsaudio.ho
jsaudio_au.o: jsaudio_au.ho
jsaudio_au.o: jsaudio_io.ho

jsaudio_dvf.o: jsaudio_dvf.c
jsaudio_dvf.o: jsaudio_io.ho
jsaudio_dvf.o: jsaudio_dvf.ho

jsaudio_io.o: jsaudio_io.c
jsaudio_io.o: jsaudio_io.ho


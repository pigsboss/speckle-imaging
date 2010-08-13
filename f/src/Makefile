CFLAGS=-L/usr/local/lib -L/usr/lib -L/home/azzuro/usr/lib -I/usr/local/include -I/usr/include -I/home/azzuro/usr/include

COPTIONS=-O3

all: average ssa deconvwnr isa si sm deconvclean smtest simobs buildsmtest rhotest

ffttest: ffttest.o libspeckle.o
	gfortran -o ffttest ffttest.o libspeckle.o smrec.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

rhotest: rhotest.o libspeckle.o
	gfortran -o rhotest rhotest.o libspeckle.o smrec.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

simobs: simobs.o libspeckle.o
	gfortran -o simobs simobs.o libspeckle.o smrec.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

buildsmtest: buildsmtest.o libspeckle.o smrec.o
	gfortran -o buildsmtest buildsmtest.o libspeckle.o smrec.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

smtest: smtest.o libspeckle.o smrec.o
	gfortran -o smtest smtest.o libspeckle.o smrec.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

sm: sm.o libspeckle.o smrec.o
	gfortran -o sm sm.o libspeckle.o smrec.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

si: si.o libspeckle.o sirec.o
	gfortran -o si si.o libspeckle.o sirec.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

isa: isa.o libspeckle.o isarec.o
	gfortran -o isa isa.o libspeckle.o isarec.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

deconvwnr: deconvwnr.o libspeckle.o
	gfortran -o deconvwnr deconvwnr.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

deconvclean: deconvclean.o libspeckle.o
	gfortran -o deconvclean deconvclean.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

ssa: ssa.o libspeckle.o ssarec.o
	gfortran -o ssa ssa.o libspeckle.o ssarec.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

average: average.o libspeckle.o
	gfortran -o average average.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

sirec.o: sirec.f
	gfortran $(COPTIONS) -c -o sirec.o sirec.f $(CFLAGS)

smrec.o: smrec.f
	gfortran $(COPTIONS) -c -o smrec.o smrec.f $(CFLAGS)

si.o: si.f
	gfortran $(COPTIONS) -c -o si.o si.f $(CFLAGS)

sm.o: sm.f
	gfortran $(COPTIONS) -c -o sm.o sm.f $(CFLAGS)

libspeckle.o: libspeckle.f
	gfortran $(COPTIONS) -c -o libspeckle.o libspeckle.f $(CFLAGS)

deconvclean.o: deconvclean.f
	gfortran $(COPTIONS) -c -o deconvclean.o deconvclean.f $(CFLAGS)

deconvwnr.o: deconvwnr.f
	gfortran $(COPTIONS) -c -o deconvwnr.o deconvwnr.f $(CFLAGS)

ssarec.o: ssarec.f
	gfortran $(COPTIONS) -c -o ssarec.o ssarec.f $(CFLAGS)

ssa.o: ssa.f
	gfortran $(COPTIONS) -c -o ssa.o ssa.f $(CFLAGS)

average.o: average.f
	gfortran $(COPTIONS) -c -o average.o average.f $(CFLAGS)

isarec.o: isarec.f
	gfortran $(COPTIONS) -c -o isarec.o isarec.f $(CFLAGS)

isa.o: isa.f
	gfortran $(COPTIONS) -c -o isa.o isa.f $(CFLAGS)

buildsmtest.o: buildsmtest.f
	gfortran $(COPTIONS) -c -o buildsmtest.o buildsmtest.f $(CFLAGS)

smtest.o: smtest.f
	gfortran $(COPTIONS) -c -o smtest.o smtest.f $(CFLAGS)

simobs.o: simobs.f
	gfortran $(COPTIONS) -c -o simobs.o simobs.f $(CFLAGS)

rhotest.o: rhotest.f
	gfortran $(COPTIONS) -c -o rhotest.o rhotest.f $(CFLAGS)

ffttest.o: ffttest.f
	gfortran $(COPTIONS) -c -o ffttest.o ffttest.f $(CFLAGS)

clean:
	rm *.o
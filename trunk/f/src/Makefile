CFLAGS=-L/usr/local/lib -L/usr/lib -L/home/azzuro/usr/lib -I/usr/local/include -I/usr/include -I/home/azzuro/usr/include

COPTIONS=-O3

all: ffc subimage bgfit bgcut average si

ffc: ffc.o libspeckle.o
	gfortran -o ffc ffc.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

subimage: subimage.o libspeckle.o
	gfortran -o subimage subimage.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

bgfit: bgfit.o libspeckle.o
	gfortran -o bgfit bgfit.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

bgcut: bgcut.o libspeckle.o
	gfortran -o bgcut bgcut.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

average: average.o libspeckle.o
	gfortran -o average average.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

si: si.o libspeckle.o
	gfortran -o si si.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

simobs: simobs.o libspeckle.o
	gfortran -o simobs simobs.o libspeckle.o -lcfitsio \
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

ffc.o: ffc.f
	gfortran $(COPTIONS) -c -o ffc.o ffc.f $(CFLAGS)

smrec.o: smrec.f
	gfortran $(COPTIONS) -c -o smrec.o smrec.f $(CFLAGS)

si.o: si.f
	gfortran $(COPTIONS) -c -o si.o si.f $(CFLAGS)

sm.o: sm.f
	gfortran $(COPTIONS) -c -o sm.o sm.f $(CFLAGS)

bgfit.o: bgfit.f
	gfortran $(COPTIONS) -c -o bgfit.o bgfit.f $(CFLAGS)

bgcut.o: bgcut.f
	gfortran $(COPTIONS) -c -o bgcut.o bgcut.f $(CFLAGS)

libspeckle.o: libspeckle.f
	gfortran $(COPTIONS) -c -o libspeckle.o libspeckle.f $(CFLAGS)

deconvclean.o: deconvclean.f
	gfortran $(COPTIONS) -c -o deconvclean.o deconvclean.f $(CFLAGS)

deconvwnr.o: deconvwnr.f
	gfortran $(COPTIONS) -c -o deconvwnr.o deconvwnr.f $(CFLAGS)

subimage.o: subimage.f
	gfortran $(COPTIONS) -c -o subimage.o subimage.f $(CFLAGS)

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

corr2.o: corr2.f
	gfortran $(COPTIONS) -c -o corr2.o corr2.f $(CFLAGS)

clean:
	rm *.o

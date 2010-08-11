all: average ssa deconvwnr isa si sm

test_sm: test_sm.o libspeckle.o smrec.o
	gfortran -o test_sm test_sm.o libspeckle.o smrec.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas

test_sm.o: test_sm.f
	gfortran -O3 -c -o test_sm.o test_sm.f -I/usr/local/include

sm: sm.o libspeckle.o smrec.o
	gfortran -o sm sm.o libspeckle.o smrec.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas

si: si.o libspeckle.o sirec.o
	gfortran -o si si.o libspeckle.o sirec.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas

isa: isa.o libspeckle.o isarec.o
	gfortran -o isa isa.o libspeckle.o isarec.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas

deconvwnr: deconvwnr.o libspeckle.o
	gfortran -o deconvwnr deconvwnr.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas

ssa: ssa.o libspeckle.o ssarec.o
	gfortran -o ssa ssa.o libspeckle.o ssarec.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas

average: average.o libspeckle.o
	gfortran -o average average.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas

sirec.o: sirec.f
	gfortran -O3 -c -o sirec.o sirec.f -I/usr/local/include

smrec.o: smrec.f
	gfortran -O3 -c -o smrec.o smrec.f -I/usr/local/include

si.o: si.f
	gfortran -O3 -c -o si.o si.f -I/usr/local/include

sm.o: sm.f
	gfortran -O3 -c -o sm.o sm.f -I/usr/local/include

libspeckle.o: libspeckle.f
	gfortran -O3 -c -o libspeckle.o libspeckle.f -I/usr/local/include

deconvwnr.o: deconvwnr.f
	gfortran -O3 -c -o deconvwnr.o deconvwnr.f -I/usr/local/include

ssarec.o: ssarec.f
	gfortran -O3 -c -o ssarec.o ssarec.f

ssa.o: ssa.f
	gfortran -O3 -c -o ssa.o ssa.f

average.o: average.f
	gfortran -O3 -c -o average.o average.f

isarec.o: isarec.f
	gfortran -O3 -c -o isarec.o isarec.f -I/usr/local/include

isa.o: isa.f
	gfortran -O3 -c -o isa.o isa.f

clean:
	rm *.o
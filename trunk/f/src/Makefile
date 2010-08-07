all: average ssa deconvwnr isa

isa: isa.o libspeckle.o isarec.o
	gfortran -o isa isa.o libspeckle.o isarec.o -lcfitsio -lacml

isarec.o: isarec.f
	gfortran -O3 -c -o isarec.o isarec.f

isa.o: isa.f
	gfortran -O3 -c -o isa.o isa.f

deconvwnr: deconvwnr.o libspeckle.o
	gfortran -o deconvwnr deconvwnr.o libspeckle.o -lcfitsio -lacml

deconvwnr.o: deconvwnr.f
	gfortran -O3 -c -o deconvwnr.o deconvwnr.f

ssa: ssa.o libspeckle.o ssarec.o
	gfortran -o ssa ssa.o libspeckle.o ssarec.o -lcfitsio -lacml

ssarec.o: ssarec.f
	gfortran -O3 -c -o ssarec.o ssarec.f

average: average.o libspeckle.o
	gfortran -o average average.o libspeckle.o -lcfitsio -lacml

libspeckle.o: libspeckle.f
	gfortran -O3 -c -o libspeckle.o libspeckle.f

ssa.o: ssa.f
	gfortran -O3 -c -o ssa.o ssa.f

average.o: average.f
	gfortran -O3 -c -o average.o average.f

clean:
	rm *.o


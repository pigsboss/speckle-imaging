all: average ssa

ssa: ssa.o libspeckle.o
	gfortran -o ssa ssa.o libspeckle.o -lcfitsio -lacml

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


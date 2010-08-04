all: average ssa

ssa: ssa.o readimage.o writeimage.o printerror.o deletefile.o
	gfortran -o ssa ssa.o readimage.o writeimage.o printerror.o deletefile.o \
	-lblas -lcfitsio

average: average_head.o average.o readimage.o writeimage.o printerror.o \
	deletefile.o
	gfortran -o average average_head.o average.o readimage.o writeimage.o \
	printerror.o deletefile.o -lblas -lcfitsio

average_head.o: average_head.f
	gfortran -O3 -c -o average_head.o average_head.f

ssa.o: ssa.f
	gfortran -O3 -c -o ssa.o ssa.f

average.o: average.f
	gfortran -O3 -c -o average.o average.f

readimage.o: readimage.f
	gfortran -O3 -c -o readimage.o readimage.f

writeimage.o: writeimage.f
	gfortran -O3 -c -o writeimage.o writeimage.f

printerror.o: printerror.f
	gfortran -O3 -c -o printerror.o printerror.f

deletefile.o: deletefile.f
	gfortran -O3 -c -o deletefile.o deletefile.f

clean:
	rm *.o

install:
	cp average ../test/
	cp ssa ../test/
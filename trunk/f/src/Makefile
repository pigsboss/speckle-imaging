CFLAGS=-L/usr/local/lib -L/usr/lib -L/home/azzuro/usr/lib -I/usr/local/include -I/usr/include -I/home/azzuro/usr/include

ALL=ffc subimage bgfit bgcut average si center ssa snr deconvwnr isa deconvclean bispostest sm spe2img img2spe fftshift ifftshift smtest kisa deconvklm

COPTIONS=-O3

all: $(ALL)

install: $(ALL)
	cp $(ALL) ../../test

deconvklm: deconvklm.o libspeckle.o
	gfortran -o deconvklm deconvklm.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

kisa: kisa.o libspeckle.o
	gfortran -o kisa kisa.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

spe2img: spe2img.o libspeckle.o
	gfortran -o spe2img spe2img.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

img2spe: img2spe.o libspeckle.o
	gfortran -o img2spe img2spe.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

fftshift: fftshift.o libspeckle.o
	gfortran -o fftshift fftshift.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

ifftshift: ifftshift.o libspeckle.o
	gfortran -o ifftshift ifftshift.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

sm: sm.o libspeckle.o
	gfortran -o sm sm.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

bispostest: bispostest.o libspeckle.o
	gfortran -o bispostest bispostest.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

snr: snr.o libspeckle.o
	gfortran -o snr snr.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

ffc: ffc.o libspeckle.o
	gfortran -o ffc ffc.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

center: center.o libspeckle.o
	gfortran -o center center.o libspeckle.o -lcfitsio \
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

ssa: ssa.o libspeckle.o
	gfortran -o ssa ssa.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

deconvwnr: deconvwnr.o libspeckle.o
	gfortran -o deconvwnr deconvwnr.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

isa: isa.o libspeckle.o
	gfortran -o isa isa.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

deconvclean: deconvclean.o libspeckle.o
	gfortran -o deconvclean deconvclean.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

smtest: smtest.o libspeckle.o
	gfortran -o smtest smtest.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

simobs: simobs.o libspeckle.o
	gfortran -o simobs simobs.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

buildsmtest: buildsmtest.o libspeckle.o smrec.o
	gfortran -o buildsmtest buildsmtest.o libspeckle.o smrec.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

img2spe.o: img2spe.f
	gfortran $(COPTIONS) -c -o img2spe.o img2spe.f $(CFLAGS)

spe2img.o: spe2img.f
	gfortran $(COPTIONS) -c -o spe2img.o spe2img.f $(CFLAGS)

ifftshift.o: ifftshift.f
	gfortran $(COPTIONS) -c -o ifftshift.o ifftshift.f $(CFLAGS)

fftshift.o: fftshift.f
	gfortran $(COPTIONS) -c -o fftshift.o fftshift.f $(CFLAGS)

kisa.o: kisa.f
	gfortran $(COPTIONS) -c -o kisa.o kisa.f $(CFLAGS)

snr.o: snr.f
	gfortran $(COPTIONS) -c -o snr.o snr.f $(CFLAGS)

ffc.o: ffc.f
	gfortran $(COPTIONS) -c -o ffc.o ffc.f $(CFLAGS)

bispostest.o: bispostest.f
	gfortran $(COPTIONS) -c -o bispostest.o bispostest.f $(CFLAGS)

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

libspeckle.o: libspeckle.f lib/math.f lib/signal.f lib/image.f lib/reconstruct.f
	gfortran $(COPTIONS) -c -o libspeckle.o libspeckle.f $(CFLAGS)

deconvklm.o: deconvklm.f
	gfortran $(COPTIONS) -c -o deconvklm.o deconvklm.f $(CFLAGS)

deconvwnr.o: deconvwnr.f
	gfortran $(COPTIONS) -c -o deconvwnr.o deconvwnr.f $(CFLAGS)

deconvclean.o: deconvclean.f
	gfortran $(COPTIONS) -c -o deconvclean.o deconvclean.f $(CFLAGS)

subimage.o: subimage.f
	gfortran $(COPTIONS) -c -o subimage.o subimage.f $(CFLAGS)

ssa.o: ssa.f
	gfortran $(COPTIONS) -c -o ssa.o ssa.f $(CFLAGS)

average.o: average.f
	gfortran $(COPTIONS) -c -o average.o average.f $(CFLAGS)

isa.o: isa.f
	gfortran $(COPTIONS) -c -o isa.o isa.f $(CFLAGS)

smtest.o: smtest.f
	gfortran $(COPTIONS) -c -o smtest.o smtest.f $(CFLAGS)

center.o: center.f
	gfortran $(COPTIONS) -c -o center.o center.f $(CFLAGS)

buildsmtest.o: buildsmtest.f
	gfortran $(COPTIONS) -c -o buildsmtest.o buildsmtest.f $(CFLAGS)

simobs.o: simobs.f
	gfortran $(COPTIONS) -c -o simobs.o simobs.f $(CFLAGS)

corr2.o: corr2.f
	gfortran $(COPTIONS) -c -o corr2.o corr2.f $(CFLAGS)

clean:
	rm *.o

CFLAGS=-L/usr/local/lib -L/usr/lib -L/home/azzuro/usr/lib -I/usr/local/include -I/usr/include -I/home/azzuro/usr/include -O3 -mtune=core2 -fPIC

ALL=ffc subimage bgfit bgcut average si center ssa snr deconvwnr isa deconvclean bispostest sm spetoimg imgtospe fftshift ifftshift smtest kisa deconvklm defringe imagesc rffm fitsinfo bmptofits simobs

FC=gfortran

all: $(ALL)

install: $(ALL)
	cp -a $(ALL) ../../test

simobs: simobs.o libspeckle.o
	$(FC) -o simobs simobs.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

bmptofits: bmptofits.o libspeckle.o
	$(FC) -o bmptofits bmptofits.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

fitsinfo: fitsinfo.o libspeckle.o
	$(FC) -o fitsinfo fitsinfo.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

rffm: rffm.o libspeckle.o
	$(FC) -o rffm rffm.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

defringe: defringe.o libspeckle.o
	$(FC) -o defringe defringe.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

imagesc: imagesc.o libspeckle.o
	$(FC) -o imagesc imagesc.o libspeckle.o -lpgplot -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

deconvklm: deconvklm.o libspeckle.o
	$(FC) -o deconvklm deconvklm.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

kisa: kisa.o libspeckle.o
	$(FC) -o kisa kisa.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

spetoimg: spetoimg.o libspeckle.o
	$(FC) -o spetoimg spetoimg.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

imgtospe: imgtospe.o libspeckle.o
	$(FC) -o imgtospe imgtospe.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

fftshift: fftshift.o libspeckle.o
	$(FC) -o fftshift fftshift.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

ifftshift: ifftshift.o libspeckle.o
	$(FC) -o ifftshift ifftshift.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

sm: sm.o libspeckle.o
	$(FC) -o sm sm.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

bispostest: bispostest.o libspeckle.o
	$(FC) -o bispostest bispostest.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

snr: snr.o libspeckle.o
	$(FC) -o snr snr.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

ffc: ffc.o libspeckle.o
	$(FC) -o ffc ffc.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

center: center.o libspeckle.o
	$(FC) -o center center.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

subimage: subimage.o libspeckle.o
	$(FC) -o subimage subimage.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

bgfit: bgfit.o libspeckle.o
	$(FC) -o bgfit bgfit.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

bgcut: bgcut.o libspeckle.o
	$(FC) -o bgcut bgcut.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

average: average.o libspeckle.o
	$(FC) -o average average.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

si: si.o libspeckle.o
	$(FC) -o si si.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

ssa: ssa.o libspeckle.o
	$(FC) -o ssa ssa.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

deconvwnr: deconvwnr.o libspeckle.o
	$(FC) -o deconvwnr deconvwnr.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

isa: isa.o libspeckle.o
	$(FC) -o isa isa.o libspeckle.o -lcfitsio -lfftw3_threads \
	-lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

deconvclean: deconvclean.o libspeckle.o
	$(FC) -o deconvclean deconvclean.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

smtest: smtest.o libspeckle.o
	$(FC) -o smtest smtest.o libspeckle.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

buildsmtest: buildsmtest.o libspeckle.o smrec.o
	$(FC) -o buildsmtest buildsmtest.o libspeckle.o smrec.o -lcfitsio \
	-lfftw3_threads -lpthread -lfftw3 -lm -llapack -lblas $(CFLAGS)

simobs.o: simobs.f
	$(FC) -c -o simobs.o simobs.f $(CFLAGS)

bmptofits.o: bmptofits.f
	$(FC) -c -o bmptofits.o bmptofits.f $(CFLAGS)

fitsinfo.o: fitsinfo.f
	$(FC) -c -o fitsinfo.o fitsinfo.f $(CFLAGS)

defringe.o: defringe.f
	$(FC) -c -o defringe.o defringe.f $(CFLAGS)

imagesc.o: imagesc.f
	$(FC) -c -o imagesc.o imagesc.f $(CFLAGS)

imgtospe.o: imgtospe.f
	$(FC) -c -o imgtospe.o imgtospe.f $(CFLAGS)

spetoimg.o: spetoimg.f
	$(FC) -c -o spetoimg.o spetoimg.f $(CFLAGS)

ifftshift.o: ifftshift.f
	$(FC) -c -o ifftshift.o ifftshift.f $(CFLAGS)

fftshift.o: fftshift.f
	$(FC) -c -o fftshift.o fftshift.f $(CFLAGS)

kisa.o: kisa.f
	$(FC) -c -o kisa.o kisa.f $(CFLAGS)

snr.o: snr.f
	$(FC) -c -o snr.o snr.f $(CFLAGS)

rffm.o: rffm.f
	$(FC) -c -o rffm.o rffm.f $(CFLAGS)

ffc.o: ffc.f
	$(FC) -c -o ffc.o ffc.f $(CFLAGS)

bispostest.o: bispostest.f
	$(FC) -c -o bispostest.o bispostest.f $(CFLAGS)

smrec.o: smrec.f
	$(FC) -c -o smrec.o smrec.f $(CFLAGS)

si.o: si.f
	$(FC) -c -o si.o si.f $(CFLAGS)

sm.o: sm.f
	$(FC) -c -o sm.o sm.f $(CFLAGS)

bgfit.o: bgfit.f
	$(FC) -c -o bgfit.o bgfit.f $(CFLAGS)

bgcut.o: bgcut.f
	$(FC) -c -o bgcut.o bgcut.f $(CFLAGS)

libspeckle.o: libspeckle.f lib/math.f lib/signal.f lib/image.f lib/reconstruct.f
	$(FC) -c -o libspeckle.o libspeckle.f $(CFLAGS)

deconvklm.o: deconvklm.f
	$(FC) -c -o deconvklm.o deconvklm.f $(CFLAGS)

deconvwnr.o: deconvwnr.f
	$(FC) -c -o deconvwnr.o deconvwnr.f $(CFLAGS)

deconvclean.o: deconvclean.f
	$(FC) -c -o deconvclean.o deconvclean.f $(CFLAGS)

subimage.o: subimage.f
	$(FC) -c -o subimage.o subimage.f $(CFLAGS)

ssa.o: ssa.f
	$(FC) -c -o ssa.o ssa.f $(CFLAGS)

average.o: average.f
	$(FC) -c -o average.o average.f $(CFLAGS)

isa.o: isa.f
	$(FC) -c -o isa.o isa.f $(CFLAGS)

smtest.o: smtest.f
	$(FC) -c -o smtest.o smtest.f $(CFLAGS)

center.o: center.f
	$(FC) -c -o center.o center.f $(CFLAGS)

buildsmtest.o: buildsmtest.f
	$(FC) -c -o buildsmtest.o buildsmtest.f $(CFLAGS)

corr2.o: corr2.f
	$(FC) -c -o corr2.o corr2.f $(CFLAGS)

clean:
	rm *.o

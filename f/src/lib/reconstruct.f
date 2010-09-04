      SUBROUTINE MEANBISP(IMGFILE,NRNG,RNG,W,ZBISP)
C  Calculate the mean bispectrum of given images.
C
C  Now only image with even NX is permitted. Otherwise BISPOS will return
C  unexpected result.
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NRNG,RNG(2,NRNG),W
      INTEGER :: STATUS
      INTEGER, PARAMETER :: UNIT=8
      DOUBLE COMPLEX, INTENT(OUT) :: ZBISP(*)
      CHARACTER(LEN=*), INTENT(IN) :: IMGFILE
C
      INTERFACE
      SUBROUTINE GETBISP(NX,NY,ZSP,Y2MAX,ZBISP)
      INTEGER, INTENT(IN) :: NX,NY,Y2MAX
      DOUBLE COMPLEX, INTENT(IN) :: ZSP(0:NX-1,0:NY-1)
      DOUBLE COMPLEX, INTENT(OUT) :: 
     &  ZBISP((Y2MAX+1)*(2*NY-Y2MAX)*NX*(NX+2)/8)
      END SUBROUTINE GETBISP
      END INTERFACE
C
C  TODO:
      STATUS=0
      RETURN
      END SUBROUTINE MEANBISP
C ******************************************************************************
      SUBROUTINE GETBISP(NX,NY,ZSP,Y2MAX,ZBISP)
C  Calculate bispectrum with given spetrum.
C
C  Now only image with even NX is permitted. Otherwise BISPOS will return
C  unexpected result.
C
      INTEGER, INTENT(IN) :: NX,NY,Y2MAX
      INTEGER :: X1,Y1,X2,Y2,K
      DOUBLE COMPLEX, INTENT(IN) :: ZSP(0:NX-1,0:NY-1)
      DOUBLE COMPLEX, INTENT(OUT) :: 
     &  ZBISP((Y2MAX+1)*(2*NY-Y2MAX)*NX*(NX+2)/8)
C
      K=1
      DO Y2=0,Y2MAX
        DO X2=0,NX-1
          DO Y1=0,NY-1-Y2
            DO X1=0,MIN(X2,NX-1-X2)
              ZBISP(K)=ZSP(X1,Y1)*ZSP(X2,Y2)*DCONJG(ZSP(X1+X2,Y1+Y2))
              K=K+1
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE GETBISP
C ******************************************************************************
      SUBROUTINE GETBISPHASE(NX,NY,DPHI,Y2MAX,DBETA)
C  Calculate bispectral phase with given spetral phase.
C
C  Now only image with even NX is permitted. Otherwise BISPOS will return
C  unexpected result.
C
      INTEGER, INTENT(IN) :: NX,NY,Y2MAX
      INTEGER :: X1,Y1,X2,Y2,K
      DOUBLE PRECISION, INTENT(IN) :: DPHI(0:NX-1,0:NY-1)
      DOUBLE PRECISION, INTENT(OUT) :: 
     &  DBETA((Y2MAX+1)*(2*NY-Y2MAX)*NX*(NX+2)/8)
C
      INTERFACE
      FUNCTION BISPOS(X1,Y1,X2,Y2,NX,NY)
      INTEGER, INTENT(IN) :: X1,Y1,X2,Y2,NX,NY
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      END INTERFACE
C
      WRITE(*,'(A,I3,A,I3)')' spectral size: ',NX,' x ',NY
      WRITE(UNIT,'(A,I3,A,I3)')' spectral size: ',NX,' x ',NY
      WRITE(*,'(A,I3)')' bispectral levels: ',Y2MAX
      WRITE(UNIT,'(A,I3)')' bispectral levels: ',Y2MAX
      K=1
      DO Y2=0,Y2MAX
        DO X2=0,NX-1
          DO Y1=0,NY-1-Y2
            DO X1=0,MIN(X2,NX-1-X2)
              DBETA(K)=DPHI(X1,Y1)+DPHI(X2,Y2)-DPHI(X1+X2,Y1+Y2)
              K=K+1
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE GETBISPHASE
C ******************************************************************************
      FUNCTION BISPOS(X1,Y1,X2,Y2,NX,NY)
C  Calculate position in bispectrum array according to positions of the
C  frequency components in the spectrum matrices.
C
C  Now only image with even NX is permitted. Otherwise the function will return
C  unexpected result.
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: X1,Y1,X2,Y2,NX,NY
      INTEGER :: BISPOS
      INTEGER :: K,L
      K=1+NX*(NX+2)*Y2*(2*NY-Y2+1)/8
      IF (X2 .LT. NX/2) THEN
        BISPOS=K+X2*(X2+1)*(NY-Y2)/2+Y1*(X2+1)+X1
        RETURN
      END IF
      IF (X2 .EQ. NX/2) THEN
        BISPOS=K+NX*(NX+2)*(NY-Y2)/8+Y1*NX/2+X1
        RETURN
      END IF
      IF (X2 .GT. NX/2) THEN
        BISPOS=K+NX*(NX+2)*(NY-Y2)/8+
     &    (3*NX-2*X2+2)*(2*X2-NX)*(NY-Y2)/8+Y1*(NX-X2)+X1
        RETURN
      END IF
      END FUNCTION BISPOS
C ******************************************************************************
      SUBROUTINE ITERATIVESHIFTADD(INFILE,NRNG,RNG,IFILE,RFILE,DSNR,
     &  NUMIT,PREFIX)
C  Iterative shift-and-add subroutine.
C
C  Arguments:
C  ==========
C
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN) :: NRNG,RNG(2,NRNG),NUMIT
      INTEGER, PARAMETER :: BUFFERSIZE=64,UNIT=8
      INTEGER :: STATUS,K,L,L1,L2,NR,LBUFFER,NBUFS,NB,NAXES(3),
     &  XM,YM,XC,YC,NFRAMES,PLANF,PLANB,INFO,NIT
      DOUBLE PRECISION, INTENT(IN) :: DSNR
      DOUBLE PRECISION :: DSUM
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:),DIMG(:,:),DREF(:,:),
     $  DEST(:,:),DCORR(:,:),DRES(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZIN(:,:),ZOUT(:,:),ZEST(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: INFILE,IFILE,RFILE,PREFIX
      CHARACTER(LEN=32) :: NUMSTR
      INTERFACE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE SHIFTADD(INFILE,NRNG,RNG,PREFIX)
      INTEGER, INTENT(IN) :: NRNG,RNG(2,NRNG)
      CHARACTER(LEN=*), INTENT(IN) :: INFILE,PREFIX
      END SUBROUTINE SHIFTADD
      SUBROUTINE DECONVWNR(NX,NY,DG,DF,DH,DSNR)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(IN) :: DG(NX,NY),DH(NX,NY),DSNR
      DOUBLE PRECISION, INTENT(OUT) :: DF(NX,NY)
      END SUBROUTINE DECONVWNR
      END INTERFACE
      STATUS=0
      OPEN(UNIT=UNIT,FILE=TRIM(PREFIX)//'_isa_summary.txt',
     &  STATUS='REPLACE',ACTION='WRITE',IOSTAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'open file failed.'
        RETURN
      END IF
      WRITE(*,'(A,I3,A)')' buffer size: ',BUFFERSIZE,'MB'
      WRITE(UNIT,'(A,I3,A)')' buffer size: ',BUFFERSIZE,'MB'
      CALL IMAGESIZE(INFILE,NAXES)
      WRITE(*,'(A,I3,A,I3)')' image size (width x height): ',
     &  NAXES(1),' x ',NAXES(2)
      WRITE(UNIT,'(A,I3,A,I3)')' image size (width x height): ',
     &  NAXES(1),' x ',NAXES(2)
      XC=INT(CEILING(0.5*REAL(NAXES(1)+1)))
      YC=INT(CEILING(0.5*REAL(NAXES(2)+1)))
      LBUFFER=NINT(DBLE(BUFFERSIZE*1024*1024)/DBLE(NAXES(1)*NAXES(2)*8))
      WRITE(*,'(A,I4,A)')' buffer length: ',LBUFFER,' frames'
      WRITE(UNIT,'(A,I4,A)')' buffer length: ',LBUFFER,' frames'
      DO NR=1,NRNG
        WRITE(*,'(A,I3,A,I5,A,I5)')' range ',NR,': from ',
     &    RNG(1,NR),' to ',RNG(2,NR)
        WRITE(UNIT,'(A,I3,A,I5,A,I5)')' range ',NR,': from ',
     &    RNG(1,NR),' to ',RNG(2,NR)
      END DO
      ALLOCATE(DBUF(NAXES(1),NAXES(2),LBUFFER),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DEST(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(IFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DEST)
      DSUM=SUM(DEST)
C     DEST=DEST*DBLE(NAXES(1)*NAXES(2))/SUM(DEST)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_isa_init.fits',
     &  (/1,1,1/),(/NAXES(1),NAXES(2),1/),DEST)
      WRITE(*,*)'initial estimate: '//TRIM(PREFIX)//'_isa_init.fits'
      WRITE(UNIT,*)'initial estimate: '//TRIM(PREFIX)//'_isa_init.fits'
      ALLOCATE(DREF(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(RFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DREF)
      DREF=DREF/SUM(DREF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_isa_psf.fits',
     &  (/1,1,1/),(/NAXES(1),NAXES(2),1/),DREF)
      WRITE(*,*)'PSF: '//TRIM(PREFIX)//'_isa_psf.fits'
      WRITE(UNIT,*)'PSF: '//TRIM(PREFIX)//'_isa_psf.fits'
      CALL DFFTW_INIT_THREADS(INFO)
      IF(INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
        RETURN
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      ALLOCATE(ZIN(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(ZOUT(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(ZEST(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DCORR(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DRES(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL DFFTW_PLAN_DFT_2D(PLANF,NAXES(1),NAXES(2),ZIN,ZOUT,-1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLANB,NAXES(1),NAXES(2),ZIN,ZOUT,1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      DO NIT=1,NUMIT
        WRITE(NUMSTR,*)NIT
        NUMSTR=TRIM(ADJUSTL(NUMSTR))
        WRITE(*,'(A,I2)')' iteration ',NIT
        CALL DECONVWNR(NAXES(1),NAXES(2),DEST,DIMG,DREF,DSNR)
        CALL WRITEIMAGE(TRIM(PREFIX)//'_isa_'//TRIM(NUMSTR)//
     &    '_core.fits',(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
        WRITE(*,*)'core '//TRIM(NUMSTR)//': '//TRIM(PREFIX)//
     &    '_isa_'//TRIM(NUMSTR)//'_core.fits'
        WRITE(UNIT,*)'core '//TRIM(NUMSTR)//': '//TRIM(PREFIX)//
     &    '_isa_'//TRIM(NUMSTR)//'_core.fits'
        ZIN=CMPLX(DIMG)
        CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
        ZEST=DCONJG(ZOUT)
        DIMG=0.0D0
        NFRAMES=0
        DO NR=1,NRNG
          NBUFS=CEILING(DBLE(RNG(2,NR)+1-RNG(1,NR))/DBLE(LBUFFER))
          DO NB=1,NBUFS
            L1=RNG(1,NR)+(NB-1)*LBUFFER
            L2=MIN(RNG(1,NR)+NB*LBUFFER-1,RNG(2,NR))
            CALL READIMAGE(INFILE,
     &       (/1,1,L1/),(/NAXES(1),NAXES(2),L2/),DBUF)
            DO L=1,L2+1-L1
              NFRAMES=NFRAMES+1
              ZIN=CMPLX(DBUF(:,:,L))
              CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
              ZIN=ZOUT*ZEST
              CALL DFFTW_EXECUTE_DFT(PLANB,ZIN,ZOUT)
              DCORR=DBLE(ZOUT)
              XM=MAXLOC(MAXVAL(DCORR,2),1)
              YM=MAXLOC(MAXVAL(DCORR,1),2)
              IF(NAXES(1)-XM .GT. XM-1)THEN
                XM=XM-1
              ELSE
                XM=XM-NAXES(1)
              END IF
              IF(NAXES(2)-YM .GT. YM-1)THEN
                YM=YM-1
              ELSE
                YM=YM-NAXES(2)
              END IF
              WRITE(*,'(A,I5,A,I2,A,I3,A,I3,A)')' estimated shifts ',
     &          NFRAMES,' in loop ',NIT,': (',XM,', ',YM,')'
              DIMG=DIMG+EOSHIFT(EOSHIFT(DBUF(:,:,L),
     &          YM,0.0D0,2),XM,0.0D0,1)
            END DO
          END DO
        END DO
        DIMG=DIMG/DBLE(NFRAMES)
        XM=MAXLOC(MAXVAL(DIMG,2),1)
        YM=MAXLOC(MAXVAL(DIMG,1),2)
        DIMG=EOSHIFT(EOSHIFT(DIMG,XM-XC,0.0D0,1),YM-YC,0.0D0,2)
        DIMG=DIMG*(DSUM/SUM(DIMG))
        CALL WRITEIMAGE(TRIM(PREFIX)//'_isa_'//TRIM(NUMSTR)//
     &    '.fits',(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
        WRITE(*,*)'output '//TRIM(NUMSTR)//': '//TRIM(PREFIX)//
     &    '_isa_'//TRIM(NUMSTR)//'.fits'
        WRITE(UNIT,*)'output '//TRIM(NUMSTR)//': '//TRIM(PREFIX)//
     &    '_isa_'//TRIM(NUMSTR)//'.fits'
        DRES=DIMG-DEST
        CALL WRITEIMAGE(TRIM(PREFIX)//'_isa_'//TRIM(NUMSTR)//
     &    '_res.fits',(/1,1,1/),(/NAXES(1),NAXES(2),1/),DRES)
        WRITE(*,*)'difference '//TRIM(NUMSTR)//': '//TRIM(PREFIX)//
     &    '_isa_'//TRIM(NUMSTR)//'_res.fits'
        WRITE(UNIT,*)'difference '//TRIM(NUMSTR)//': '//TRIM(PREFIX)//
     &    '_isa_'//TRIM(NUMSTR)//'_res.fits'
        WRITE(*,'(A,ES10.3)')' total difference '//NUMSTR//': ',
     &    DSQRT(SUM(DRES*DRES)/DBLE(NAXES(1)*NAXES(2)))
        WRITE(UNIT,'(A,ES10.3)')' total difference '//NUMSTR//': ',
     &    DSQRT(SUM(DRES*DRES)/DBLE(NAXES(1)*NAXES(2)))
        DEST=DIMG
      END DO
      CALL DFFTW_DESTROY_PLAN(PLANF)
      CALL DFFTW_DESTROY_PLAN(PLANB)
      CLOSE(UNIT)
      DEALLOCATE(DBUF)
      DEALLOCATE(DIMG)
      DEALLOCATE(DREF)
      DEALLOCATE(DEST)
      DEALLOCATE(DCORR)
      DEALLOCATE(ZIN)
      DEALLOCATE(DRES)
      DEALLOCATE(ZOUT)
      DEALLOCATE(ZEST)
      RETURN
      END SUBROUTINE ITERATIVESHIFTADD
C ******************************************************************************
      SUBROUTINE SHIFTADD(INFILE,NRNG,RNG,PREFIX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NRNG,RNG(2,NRNG)
      INTEGER, PARAMETER :: BUFFERSIZE=64
      INTEGER :: STATUS,K,L,L1,L2,NR,LBUFFER,NBUFS,NB,NAXES(3),
     &  XM,YM,XC,YC,NFRAMES
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:),DIMG(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: INFILE,PREFIX
      INTERFACE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      END INTERFACE
      STATUS=0
      WRITE(*,'(A,I3,A)')' buffer size: ',BUFFERSIZE,'MB'
      CALL IMAGESIZE(INFILE,NAXES)
      XC=INT(CEILING(0.5*REAL(1+NAXES(1))))
      YC=INT(CEILING(0.5*REAL(1+NAXES(2))))
      WRITE(*,'(A,I3,A,I3)')' image size (width x height): ',
     &  NAXES(1),' x ',NAXES(2)
      LBUFFER=NINT(DBLE(BUFFERSIZE*1024*1024)/DBLE(NAXES(1)*NAXES(2)*8))
      WRITE(*,'(A,I4,A)')' buffer length: ',LBUFFER,' frames'
      DO NR=1,NRNG
        WRITE(*,'(A,I3,A,I5,A,I5)')' range ',NR,': from ',
     &    RNG(1,NR),' to ',RNG(2,NR)
      END DO
      ALLOCATE(DBUF(NAXES(1),NAXES(2),LBUFFER),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      DIMG=0.0D0
      NFRAMES=0
      DO NR=1,NRNG
        NBUFS=CEILING(DBLE(RNG(2,NR)+1-RNG(1,NR))/DBLE(LBUFFER))
        DO NB=1,NBUFS
          L1=RNG(1,NR)+(NB-1)*LBUFFER
          L2=MIN(RNG(1,NR)+NB*LBUFFER-1,RNG(2,NR))
          CALL READIMAGE(INFILE,
     &     (/1,1,L1/),(/NAXES(1),NAXES(2),L2/),DBUF)
          DO L=1,L2+1-L1
            NFRAMES=NFRAMES+1
            XM=MAXLOC(MAXVAL(DBUF(:,:,L),2),1)
            YM=MAXLOC(MAXVAL(DBUF(:,:,L),1),2)
            WRITE(*,'(A,I5,A,I3,A,I3,A)')' maximum location ',
     &        NFRAMES,': (',XM,', ',YM,')'
            DIMG=DIMG+EOSHIFT(EOSHIFT(DBUF(:,:,L),
     &        YM-YC,0.0D0,2),XM-XC,0.0D0,1)
          END DO
        END DO
      END DO
      DIMG=DIMG/DBLE(NFRAMES)
      DEALLOCATE(DBUF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_ssa.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DIMG)
      PRINT *,'output: '//TRIM(PREFIX)//'_ssa.fits'
      DEALLOCATE(DIMG)
      RETURN
      END SUBROUTINE SHIFTADD
C ******************************************************************************
      SUBROUTINE SPECKLEINTERFEROMETRY(TFILE,TRNG,RFILE,RRNG,PX,PY,
     &  PREFIX,BUFFERSIZE)
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER, INTENT(IN) :: PX,PY,TRNG(2),RRNG(2)
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      INTEGER :: X,Y,NAXES(3)
      DOUBLE PRECISION, DIMENSION(PX,PY) :: DTPSD,DRPSD,DPSD,DACF
      DOUBLE COMPLEX :: ZIN(PX,PY),ZOUT(PX,PY)
      CHARACTER*(*), INTENT(IN) :: TFILE,RFILE,PREFIX
      INTERFACE
      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      INTEGER, INTENT(OUT) :: NAXES(3)
      END SUBROUTINE IMAGESIZE
      SUBROUTINE GETPSD(FILENAME,FPIXELS,LPIXELS,DX,DY,DPSD,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3),DX,DY
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      DOUBLE PRECISION, INTENT(OUT) :: DPSD(DX,DY)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE GETPSD
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: ARRAY(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
      SUBROUTINE DFFTSHIFT(W,H,DX)
      INTEGER, INTENT(IN) :: W,H
      DOUBLE PRECISION, INTENT(INOUT) :: DX(W,H)
      END SUBROUTINE DFFTSHIFT
      SUBROUTINE DIFFTSHIFT(W,H,DX)
      INTEGER, INTENT(IN) :: W,H
      DOUBLE PRECISION, INTENT(INOUT) :: DX(W,H)
      END SUBROUTINE DIFFTSHIFT
      END INTERFACE
C  Statements:
C  ===========
      CALL IMAGESIZE(TFILE,NAXES)
      IF(TRNG(2).LT.NAXES(3))THEN
        PRINT *,'request is out of range.'
        RETURN
      END IF
      IF(PRESENT(BUFFERSIZE))THEN
        CALL GETPSD(TFILE,(/1,1,TRNG(1)/),(/NAXES(1),NAXES(2),TRNG(2)/),
     &    PX,PY,DTPSD,BUFFERSIZE)
      ELSE
        CALL GETPSD(TFILE,(/1,1,TRNG(1)/),(/NAXES(1),NAXES(2),TRNG(2)/),
     &    PX,PY,DTPSD)
      END IF
      CALL DFFTSHIFT(NAXES(1),NAXES(2),DTPSD)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_tgabs.fits',
     &  (/1,1,1/),(/PX,PY,1/),DSQRT(DTPSD))
      PRINT *,'spectral amplitude of target: '//
     &  TRIM(PREFIX)//'_tgabs.fits'
      CALL IMAGESIZE(RFILE,NAXES)
      IF(RRNG(2).LT.NAXES(3))THEN
        PRINT *,'request is out of range.'
        RETURN
      END IF
      IF(PRESENT(BUFFERSIZE))THEN
        CALL GETPSD(RFILE,(/1,1,RRNG(1)/),(/NAXES(1),NAXES(2),RRNG(2)/),
     &    PX,PY,DRPSD,BUFFERSIZE)
      ELSE
        CALL GETPSD(RFILE,(/1,1,RRNG(1)/),(/NAXES(1),NAXES(2),RRNG(2)/),
     &    PX,PY,DRPSD)
      END IF
      CALL DFFTSHIFT(NAXES(1),NAXES(2),DRPSD)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_rfabs.fits',
     &  (/1,1,1/),(/PX,PY,1/),DSQRT(DRPSD))
      PRINT *,'spectral amplitude of reference: '//
     &  TRIM(PREFIX)//'_rfabs.fits'
      DO X=1,PX
        DO Y=1,PY
          IF (DABS(DRPSD(X,Y)) .GT. 0)THEN
            DPSD(X,Y)=DTPSD(X,Y)/DRPSD(X,Y)
          ELSE
            DPSD(X,Y)=DBLE(0)
          END IF
        END DO
      END DO
      CALL WRITEIMAGE(TRIM(PREFIX)//'_dmpsd.fits',
     &  (/1,1,1/),(/PX,PY,1/),DPSD)
      PRINT *,'demodulated power spectral density: '//
     &  TRIM(PREFIX)//'_dmpsd.fits'
      CALL DFFTW_PLAN_DFT_2D(PLAN,PX,PY,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      CALL DIFFTSHIFT(NAXES(1),NAXES(2),DPSD)
      ZIN=CMPLX(DPSD)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DACF=DREAL(ZOUT/DBLE(NAXES(1))/DBLE(NAXES(2)))
      CALL DFFTSHIFT(PX,PY,DACF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_acf.fits',
     &  (/1,1,1/),(/PX,PY,1/),DACF)
      PRINT *,'reconstructed auto-correlation function: '//
     &  TRIM(PREFIX)//'_acf.fits'
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE SPECKLEINTERFEROMETRY
C ******************************************************************************
      SUBROUTINE GETPSD(FILENAME,FPIXELS,LPIXELS,DX,DY,DPSD,BUFFERSIZE)
C  Purpose:
C  ========
C  Get the mean power spectral density of all frames in the given FITS file.
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN):: FPIXELS(3),LPIXELS(3),DX,DY
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      INTEGER :: W,H,INFO,K,NPIXELS,NFRAMES,LBUF,NBUF,L,L1,L2
      INTEGER*8 :: PLAN
      DOUBLE PRECISION, INTENT(OUT) :: DPSD(DX,DY)
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:)
      DOUBLE COMPLEX, DIMENSION(DX,DY) :: ZIN,ZOUT
      CHARACTER*(*), INTENT(IN) :: FILENAME
      INTERFACE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      END INTERFACE
C  Statements:
C  ===========
      W=LPIXELS(1)-FPIXELS(1)+1
      H=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=W*H
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      IF(PRESENT(BUFFERSIZE))THEN
        LBUF=INT(FLOOR(DBLE(BUFFERSIZE*1024*1024)/DBLE(NPIXELS*8)))
      ELSE
        LBUF=INT(FLOOR(DBLE(1024*1024)/DBLE(NPIXELS)))
      END IF
      NBUF=INT(CEILING(DBLE(NFRAMES)/DBLE(LBUF)))
      ALLOCATE(DBUF(W,H,LBUF))
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
        RETURN
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_IMPORT_SYSTEM_WISDOM(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_IMPORT_SYSTEM_WISDOM failed.'
      END IF
C     PRINT *,'Start planning.'
      CALL DFFTW_PLAN_DFT_2D(PLAN,DX,DY,ZIN,ZOUT,-1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
C     PRINT *,'Finished planning.'
      DPSD=DBLE(0)
      DO K=1,NBUF
        L1=(K-1)*LBUF+1
        L2=MIN(NFRAMES,K*LBUF)
        CALL READIMAGE(FILENAME,(/FPIXELS(1),FPIXELS(2),L1/),
     &    (/LPIXELS(1),LPIXELS(2),L2/),DBUF)
        DO L=1,L2+1-L1
          ZIN=CMPLX(DBLE(0))
          ZIN(1:W,1:H)=DBUF(:,:,L)*DBLE(NPIXELS)/SUM(DBUF(:,:,L))
          CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
          DPSD=DPSD+DREAL(ZOUT*CONJG(ZOUT))
        END DO
      END DO
      DPSD=DPSD/DBLE(NFRAMES)
      DEALLOCATE(DBUF)
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE GETPSD
C ******************************************************************************
      SUBROUTINE DECONVCLEAN(NX,NY,DG,DF,DH,DBETA,NUMIT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX,NY,NUMIT
      INTEGER :: X,Y,XC,YC,K,L
      DOUBLE PRECISION, INTENT(IN) :: DG(NX,NY),DH(NX,NY),DBETA
      DOUBLE PRECISION, INTENT(OUT) :: DF(NX,NY)
      DOUBLE PRECISION :: DRES(NX,NY),DPSF(NX,NY),DTMP(NX,NY),DSTD
      INTERFACE
      END INTERFACE
      XC=INT(CEILING(0.5*DBLE(NX+1)))
      YC=INT(CEILING(0.5*DBLE(NY+1)))
      X=MAXLOC(MAXVAL(DH,2),1)
      Y=MAXLOC(MAXVAL(DH,1),2)
      DPSF=EOSHIFT(EOSHIFT(DH,Y-YC,0.0D0,2),X-XC,0.0D0,1)
      DPSF=DPSF/SUM(DPSF)
      DRES=DG
      L=1
      DO K=1,NUMIT
        L=L+1
        X=MAXLOC(MAXVAL(DRES,2),1)
        Y=MAXLOC(MAXVAL(DRES,1),2)
        DF(X,Y)=DF(X,Y)+DRES(X,Y)*DBETA
        DTMP=EOSHIFT(EOSHIFT(DPSF,YC-Y,0.0D0,2),XC-X,0.0D0,1)
        DTMP=DTMP/SUM(DTMP)
        DRES=DRES-DRES(X,Y)*DBETA*DTMP
        IF (SUM(DRES) .LE. 0.0D0)THEN
          EXIT
        END IF
        IF (L .GT. 1000)THEN
          L=1
          PRINT *,'Loop: ',K
          DSTD=DSQRT(SUM((DRES-SUM(DRES)/DBLE(NX*NY))*
     &      (DRES-SUM(DRES)/DBLE(NX*NY)))/DBLE(NX*NY))
          WRITE (*,'(A,ES9.2,A,ES9.2)') ' deviation of residual: ',DSTD,
     &      ', maximum of PSF: ',MAXVAL(DPSF)
          WRITE (*,'(A,ES10.3,A,ES10.3,A,ES10.3)')
     &     ' Residual sum:',SUM(DRES),
     &      ', CLEAN sum:',SUM(DF),', Total sum:',SUM(DF)+SUM(DRES)
          WRITE (*,'(A,I3,A,I3,A,ES10.3,ES10.3)')
     &      ' Maxima: (',X,',',Y,'),',MAXVAL(DRES),DRES(Y,X)
          IF(DSTD .LE. MAXVAL(DPSF))THEN
            PRINT *,'CLEAN deep enough, stop.'
            EXIT
          END IF
        END IF
      END DO
      CALL WRITEIMAGE('CLEAN_RESIDUAL.FITS',(/1,1,1/),(/NX,NY,1/),DRES)
      RETURN
      END SUBROUTINE DECONVCLEAN
C ******************************************************************************
      SUBROUTINE DECONVWNR(NX,NY,DG,DF,DH,DSNR)
C  Wiener deconvolution subroutine.
C
C  Purpose:
C  ========
C  G = conv(F, H). This routine returns F.
C
C  Arguments:
C  ==========
C  NX - Number of rows of G.
C  NY - Number of columns of G.
C  DG - Input.
C  DF - Output. size(F) = size(G).
C  DH - Input. size(PSF) = size (G)
C  DSNR - Signal-to-noise ratio.
C 
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER, INTENT(IN) :: NX,NY
      INTEGER :: INFO,X,Y
      DOUBLE PRECISION, INTENT(IN) :: DG(NX,NY),DH(NX,NY),DSNR
      DOUBLE PRECISION, INTENT(OUT) :: DF(NX,NY)
      DOUBLE PRECISION :: DTMP
      DOUBLE COMPLEX :: ZG(NX,NY),ZH(NX,NY)
      DOUBLE COMPLEX :: ZIN(NX,NY),ZOUT(NX,NY),ZDECONV(NX,NY)
      INTERFACE
      SUBROUTINE ZFFTSHIFT(NX,NY,ZX)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE COMPLEX, INTENT(INOUT) :: ZX(NX,NY)
      END SUBROUTINE ZFFTSHIFT
      END INTERFACE
      CALL DFFTW_PLAN_DFT_2D(PLAN,NX,NY,ZIN,ZOUT,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=CMPLX(DG)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZG=ZOUT
      ZIN=CMPLX(DH)
      CALL ZFFTSHIFT(NX,NY,ZIN)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZH=ZOUT
      ZDECONV=DCONJG(ZH)/(ZH*DCONJG(ZH)+CMPLX(DBLE(1)/DSNR))
      CALL DFFTW_PLAN_DFT_2D(PLAN,NX,NY,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZG*ZDECONV
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DF=DBLE(ZOUT)/DBLE(NX*NY)
      DTMP=SUM(DF)
      DO X=1,NX
        DO Y=1,NY
          IF(DF(X,Y) .LT. 0.0D0)THEN
            DF(X,Y)=0.0D0
          END IF
        END DO
      END DO
      DF=DF*DTMP/SUM(DF)
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE DECONVWNR
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************

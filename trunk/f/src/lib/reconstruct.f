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
      XC=INT(FLOOR(0.5*REAL(1+NAXES(1))))
      YC=INT(FLOOR(0.5*REAL(1+NAXES(2))))
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
      SUBROUTINE DECONVCLEAN(M,N,DG,DF,DH,DBETA,MNUMIT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: X,Y,XC,YC,K,L,MNUMIT
      DOUBLE PRECISION, INTENT(IN) :: DG(M,N),DH(M,N),DBETA
      DOUBLE PRECISION, INTENT(OUT) :: DF(M,N)
      DOUBLE PRECISION :: DRES(M,N),DPSF(M,N),DSNR,DR
      XC=INT(CEILING(0.5*DBLE(N+1)))
      YC=INT(CEILING(0.5*DBLE(M+1)))
      DR=DBLE(0.5)*DBLE(MIN(M,N))
      X=MAXLOC(MAXVAL(DH,1),2)
      Y=MAXLOC(MAXVAL(DH,2),1)
      DPSF=EOSHIFT(EOSHIFT(DH,Y-YC,DBLE(0),1),X-XC,DBLE(0),2)
      DPSF=DPSF/SUM(DPSF)
      DRES=DG
      L=1
      DO K=1,MNUMIT
        L=L+1
        X=MAXLOC(MAXVAL(DRES,1),2)
        Y=MAXLOC(MAXVAL(DRES,2),1)
        DF(Y,X)=DF(Y,X)+DRES(Y,X)*DBETA
        DRES=DRES-DRES(Y,X)*DBETA*EOSHIFT(EOSHIFT(DPSF,YC-Y,DBLE(0),1),
     &    XC-X,DBLE(0),2)
        IF (SUM(DRES) .LE. 0)THEN
          EXIT
        END IF
        IF (L .GT. 100)THEN
          L=1
          PRINT *,'Loop: ',K
          CALL GETSNR(M,N,DR,DSNR,DRES,'P4')
          WRITE (*,'(A,ES10.3)') ' Residual SNR:',DSNR
          WRITE (*,'(A,ES10.3,A,ES10.3,A,ES10.3)')
     &     ' Residual sum:',SUM(DRES),
     &      ', CLEAN sum:',SUM(DF),', Total sum:',SUM(DF)+SUM(DRES)
          WRITE (*,'(A,I3,A,I3,A,ES10.3,ES10.3)')
     &      ' Maxima: (',X,',',Y,'),',MAXVAL(DRES),DRES(Y,X)
        END IF
      END DO
      CALL WRITEIMAGE('CLEAN_RESIDUAL.FITS',(/1,1,1/),(/M,N,1/),DRES)
      RETURN
      END SUBROUTINE DECONVCLEAN
C ******************************************************************************
      SUBROUTINE DECONVWNR(M,N,DG,DF,DH,DSNR)
C  Wiener deconvolution subroutine.
C
C  Purpose:
C  ========
C  G = conv(F, H). This routine returns F.
C
C  Arguments:
C  ==========
C  M - Number of rows of G.
C  N - Number of columns of G.
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
      INTEGER :: M,N,INFO
      DOUBLE PRECISION :: DG(M,N),DF(M,N),DH(M,N)
      DOUBLE PRECISION :: DSNR
      DOUBLE COMPLEX :: ZG(M,N),ZH(M,N)
      DOUBLE COMPLEX :: ZIN(M,N),ZOUT(M,N),ZDECONV(M,N)
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=CMPLX(DG)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZG=ZOUT
      ZIN=CMPLX(DH)
      CALL ZFFTSHIFT(M,N,ZIN)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZH=ZOUT
      ZDECONV=CONJG(ZH)/(ZH*CONJG(ZH)+CMPLX(DBLE(1)/DSNR))
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZG*ZDECONV
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DF=DBLE(ZOUT)
      CALL DFFTW_DESTROY_PLAN(PLAN)
C     CALL DFFTW_CLEANUP_THREADS()
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

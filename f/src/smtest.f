      PROGRAM SMTEST
      IMPLICIT NONE
      INTEGER :: STATUS,PLAN,NAXES(3),LBISP,Y2MAX,NARGS,K
      INTEGER, PARAMETER :: DEFAULTY2MAX=3,UNIT=8,NUMFRMS=100
      DOUBLE PRECISION, ALLOCATABLE :: DIMG(:,:),DRHO(:,:),DPHI(:,:),
     &  DBISP(:),DTMP(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZSP(:,:),ZBISP(:)
      CHARACTER(LEN=256) :: FILENAME,PREFIX,ARG
C  interface:
C  ==========
      INTERFACE
      SUBROUTINE MEANBISP(IMGFILE,NRNG,RNG,Y2MAX,ZBISP)
      INTEGER, INTENT(IN) :: NRNG,RNG(2,NRNG),Y2MAX
      DOUBLE COMPLEX, INTENT(OUT) :: ZBISP(*)
      CHARACTER(LEN=*), INTENT(IN) :: IMGFILE
      END SUBROUTINE MEANBISP
      SUBROUTINE GETARGUMENT(NX,NY,ZX,DPHI)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(OUT) :: DPHI(NX,NY)
      DOUBLE COMPLEX, INTENT(IN) :: ZX(NX,NY)
      END SUBROUTINE GETARGUMENT
      SUBROUTINE RECURSPHASE(NX,NY,Y2MAX,DBETA,DPHI)
      INTEGER, INTENT(IN) :: NX,NY,Y2MAX
      DOUBLE PRECISION, INTENT(OUT) :: DPHI(0:NX-1,0:NY-1)
      DOUBLE PRECISION, INTENT(IN) :: 
     &  DBETA((Y2MAX+1)*(2*NY-Y2MAX)*NX*(NX+2)/8)
      END SUBROUTINE RECURSPHASE
      SUBROUTINE ADDBISP(NX,NY,ZSP,Y2MAX,ZBISP)
      INTEGER, INTENT(IN) :: NX,NY,Y2MAX
      DOUBLE COMPLEX, INTENT(IN) :: ZSP(0:NX-1,0:NY-1)
      DOUBLE COMPLEX, INTENT(INOUT) :: 
     &  ZBISP((Y2MAX+1)*(2*NY-Y2MAX)*NX*(NX+2)/8)
      END SUBROUTINE ADDBISP
      FUNCTION BISPOS(X1,Y1,X2,Y2,NX,NY)
      INTEGER, INTENT(IN) :: X1,Y1,X2,Y2,NX,NY
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      SUBROUTINE SPECTRUMTOIMAGE(NX,NY,ZSP,DIMG)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(NX,NY)
      DOUBLE COMPLEX, INTENT(IN) :: ZSP(NX,NY)
      END SUBROUTINE SPECTRUMTOIMAGE
      SUBROUTINE IMAGETOSPECTRUM(NX,NY,DIMG,ZSP)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(IN) :: DIMG(NX,NY)
      DOUBLE COMPLEX, INTENT(OUT) :: ZSP(NX,NY)
      END SUBROUTINE IMAGETOSPECTRUM
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
C  end interface
C  =============
      STATUS=0
      CALL GET_COMMAND_ARGUMENT(1,FILENAME)
      PREFIX='smtest'
      OPEN(UNIT=UNIT,FILE=TRIM(PREFIX)//'.log',STATUS='REPLACE',
     &  ACTION='WRITE',IOSTAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'open file failed.'
        RETURN
      END IF
      Y2MAX=DEFAULTY2MAX
C
      NARGS=COMMAND_ARGUMENT_COUNT()
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-level=') .GT. 0)THEN
          READ(ARG(INDEX(ARG,'-level=')+7:),*)Y2MAX
        END IF
      END DO
C
      CALL IMAGESIZE(FILENAME,NAXES)
      LBISP=(Y2MAX+1)*(2*NAXES(1)-Y2MAX)*NAXES(1)*(NAXES(1)+2)/8
      WRITE(*,*)'length of bispectral array: ',LBISP
      WRITE(UNIT,*)'length of bispectral array: ',LBISP
      
C
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        RETURN
      END IF
      ALLOCATE(DTMP(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        RETURN
      END IF
      ALLOCATE(ZSP(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        RETURN
      END IF
      ALLOCATE(DRHO(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        RETURN
      END IF
      ALLOCATE(DPHI(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        RETURN
      END IF
      ALLOCATE(DBISP(LBISP),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        RETURN
      END IF
      ALLOCATE(ZBISP(LBISP),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        RETURN
      END IF
C
      CALL READIMAGE(FILENAME,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      CALL IMAGETOSPECTRUM(NAXES(1),NAXES(2),DIMG,ZSP)
C
      CALL WRITEIMAGE(TRIM(PREFIX)//'_i_real.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DREAL(ZSP))
      WRITE(*,*)'real part of spectrum: '//
     &  TRIM(PREFIX)//'_i_real.fits'
      WRITE(UNIT,*)'real part of spectrum: '//
     &  TRIM(PREFIX)//'_i_real.fits'
      CALL WRITEIMAGE(TRIM(PREFIX)//'_i_imag.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DIMAG(ZSP))
      WRITE(*,*)'imaginary part of spectrum: '//
     &  TRIM(PREFIX)//'_i_imag.fits'
      WRITE(UNIT,*)'imaginary part of spectrum: '//
     &  TRIM(PREFIX)//'_i_imag.fits'
C
      DRHO=ZABS(ZSP)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_i_mod.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DRHO)
      WRITE(*,*)'input spectral modulus: '//
     &  TRIM(PREFIX)//'_i_mod.fits'
      WRITE(UNIT,*)'input spectral modulus: '//
     &  TRIM(PREFIX)//'_i_mod.fits'
C
      CALL GETARGUMENT(NAXES(1),NAXES(2),ZSP,DPHI)
      DTMP=DPHI
      DPHI=DATAN2(DIMAG(ZSP),DREAL(ZSP))
      WRITE(*,'(A,ES8.2)')' norm of phase difference: ',
     &  DSQRT(SUM(DTMP*DTMP)/DBLE(NAXES(1)*NAXES(2)))
      WRITE(UNIT,'(A,ES8.2)')' norm of phase difference: ',
     &  DSQRT(SUM(DTMP*DTMP)/DBLE(NAXES(1)*NAXES(2)))
C
      CALL WRITEIMAGE(TRIM(PREFIX)//'_i_pha.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DPHI)
      WRITE(*,*)'input spectral phase: '//
     &  TRIM(PREFIX)//'_i_pha.fits'
      WRITE(UNIT,*)'input spectral phase: '//
     &  TRIM(PREFIX)//'_i_pha.fits'
C
      ZSP=DRHO*DCMPLX(DCOS(DPHI),DSIN(DPHI))
      DTMP=DIMG
      CALL SPECTRUMTOIMAGE(NAXES(1),NAXES(2),ZSP,DIMG)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_i_img.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DIMG)
      WRITE(*,*)'input image recovered from spectrum: '//
     &  TRIM(PREFIX)//'_i_img.fits'
      WRITE(UNIT,*)'input image recovered from spectrum: '//
     &  TRIM(PREFIX)//'_i_img.fits'
C
      DTMP=DTMP-DIMG
      CALL WRITEIMAGE(TRIM(PREFIX)//'_i_img_diff.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DTMP)
      WRITE(*,*)'difference introduced by fourier transform: '//
     &  TRIM(PREFIX)//'_i_img_diff.fits'
      WRITE(UNIT,*)'difference introduced by fourier transform: '//
     &  TRIM(PREFIX)//'_i_img_diff.fits'
      WRITE(*,'(A,ES8.2)')' norm of difference: ',
     &  DSQRT(SUM(DTMP*DTMP)/DBLE(NAXES(1)*NAXES(2)))
      WRITE(UNIT,'(A,ES8.2)')' norm of difference: ',
     &  DSQRT(SUM(DTMP*DTMP)/DBLE(NAXES(1)*NAXES(2)))
      DIMG=DTMP
C
      DTMP=DPHI
      DPHI=DSIN(DPHI)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_i_sin.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DPHI)
      WRITE(*,*)'sine of input spectral phase: '//
     &  TRIM(PREFIX)//'_i_sin.fits'
      WRITE(UNIT,*)'sine of input spectral phase: '//
     &  TRIM(PREFIX)//'_i_sin.fits'
C
      ZBISP=DCMPLX(0.0D0)
C     CALL ADDBISP(NAXES(1),NAXES(2),ZSP,Y2MAX,ZBISP)
      DO K=1,NUMFRMS
        CALL ADDBISP(NAXES(1),NAXES(2),ZSP,Y2MAX,ZBISP)
      END DO
      ZBISP=ZBISP/DBLE(NUMFRMS)
C     CALL MEANBISP(FILENAME,1,(/1,1/),Y2MAX,ZBISP)
C     CALL GETARGUMENT(LBISP,1,ZBISP,DBISP)
      DBISP=DATAN2(DIMAG(ZBISP),DREAL(ZBISP))
      CALL RECURSPHASE(NAXES(1),NAXES(2),Y2MAX,DBISP,DPHI)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_o_pha.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DPHI)
      WRITE(*,*)'phase recured from bispectrum: '//
     &  TRIM(PREFIX)//'_o_pha.fits'
      WRITE(UNIT,*)'phase recured from bispectrum: '//
     &  TRIM(PREFIX)//'_o_pha.fits'
C
      ZSP=DRHO*DCMPLX(DCOS(DPHI),DSIN(DPHI))
      DTMP=DIMG
      CALL SPECTRUMTOIMAGE(NAXES(1),NAXES(2),ZSP,DIMG)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_o_img.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DIMG)
      WRITE(*,*)'input image recovered from bispectrum: '//
     &  TRIM(PREFIX)//'_o_img.fits'
      WRITE(UNIT,*)'input image recovered from bispectrum: '//
     &  TRIM(PREFIX)//'_o_img.fits'
C
      DIMG=DTMP-DIMG
      CALL WRITEIMAGE(TRIM(PREFIX)//'_o_img_diff.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DIMG)
      WRITE(*,*)'difference introduced by recursion: '//
     &  TRIM(PREFIX)//'_o_img_diff.fits'
      WRITE(UNIT,*)'difference introduced by recursion: '//
     &  TRIM(PREFIX)//'_o_img_diff.fits'
      WRITE(*,'(A,ES8.2)')' norm of difference: ',
     &  DSQRT(SUM(DIMG*DIMG)/DBLE(NAXES(1)*NAXES(2)))
      WRITE(UNIT,'(A,ES8.2)')' norm of difference: ',
     &  DSQRT(SUM(DIMG*DIMG)/DBLE(NAXES(1)*NAXES(2)))
C
      DEALLOCATE(DIMG)
      DEALLOCATE(DTMP)
      DEALLOCATE(ZSP)
      DEALLOCATE(DRHO)
      DEALLOCATE(DPHI)
      DEALLOCATE(DBISP)
      DEALLOCATE(ZBISP)
      CLOSE(UNIT)
      STOP
      END PROGRAM SMTEST

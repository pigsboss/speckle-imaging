      SUBROUTINE IMFILTER(NX,NY,DF,DH,DG)
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(OUT) :: DF(NX,NY),DH(NX,NY),DG(NX,NY)
      INTEGER :: PLANF,PLANB
      DOUBLE COMPLEX :: ZIN(NX,NY),ZOUT(NX,NY),ZH(NX,NY)
      INTERFACE
      SUBROUTINE DFFTSHIFT(NX,NY,DIMG)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(INOUT) :: DIMG(NX,NY)
      END SUBROUTINE DFFTSHIFT
      SUBROUTINE DIFFTSHIFT(NX,NY,DIMG)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(INOUT) :: DIMG(NX,NY)
      END SUBROUTINE DIFFTSHIFT
      SUBROUTINE ZFFTSHIFT(NX,NY,ZSP)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE COMPLEX, INTENT(INOUT) :: ZSP(NX,NY)
      END SUBROUTINE ZFFTSHIFT
      SUBROUTINE ZIFFTSHIFT(NX,NY,ZSP)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE COMPLEX, INTENT(INOUT) :: ZSP(NX,NY)
      END SUBROUTINE ZIFFTSHIFT
      END INTERFACE
      CALL DFFTW_PLAN_DFT_2D(PLANF,NX,NY,ZIN,ZOUT,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLANB,NX,NY,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=DCMPLX(DH)
      CALL ZIFFTSHIFT(NX,NY,ZIN)
      CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZH)
      ZIN=DCMPLX(DF)
      CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
      ZIN=ZH*ZOUT/DBLE(NX)/DBLE(NY)
      CALL DFFTW_EXECUTE_DFT(PLANB,ZIN,ZOUT)
      DG=DREAL(ZOUT)/DSQRT(DBLE(NX)*DBLE(NY))
      CALL DFFTW_DESTROY_PLAN(PLANF)
      CALL DFFTW_DESTROY_PLAN(PLANB)
      RETURN
      END SUBROUTINE IMFILTER
C ******************************************************************************
      SUBROUTINE FLATFIELDCORRECT(TARGFILE,FLATFILE,DARKFILE,
     &  FPIXELS,LPIXELS,CRCTFILE)
C  Flat-field correction subroutine.
C
C  Usage:
C  =====
C  Input direct observation results, flat field file, as well as dark current
C  file. Return flat-field corrected observation results as output.
C
C                 I_original - I_dark
C  I_corrected = ---------------------
C                   I_flat - I_dark
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: STATUS,UNIT,RWMODE,BLOCKSIZE,GROUP
      INTEGER :: K,M,N,NPIXELS,NFRAMES,NAXES(3),NFOUND
      DOUBLE PRECISION, ALLOCATABLE :: DFLAT(:),DDARK(:),DTARG(:),
     &  DCRCT(:)
      CHARACTER*(*), INTENT(IN) :: TARGFILE,FLATFILE,DARKFILE,CRCTFILE
      CHARACTER(LEN=256) :: KEYVAL,COMMENT
      INTERFACE
      SUBROUTINE AVERAGE(AVGFILE,FPIXELS,LPIXELS,DAVG,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER, INTENT(INOUT), OPTIONAL :: BUFFERSIZE
      DOUBLE PRECISION :: DAVG(LPIXELS(1)+1-FPIXELS(1),
     &  LPIXELS(2)+1-FPIXELS(2))
      CHARACTER*(*), INTENT(IN) :: AVGFILE
      END SUBROUTINE AVERAGE
      END INTERFACE
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      RWMODE=0
      CALL FTOPEN(UNIT,TARGFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(TARGFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      M=NAXES(1)
      N=NAXES(2)
      CALL FTCLOS(UNIT, STATUS)
      CALL FTOPEN(UNIT,FLATFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(FLATFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      IF ((NAXES(1).NE.M) .OR. (NAXES(2).NE.N))THEN
        PRINT *,TRIM(TARGFILE)//' and '//TRIM(FLATFILE)//
     &    ' differ in size.'
        RETURN
      END IF
      CALL FTCLOS(UNIT, STATUS)
      ALLOCATE(DFLAT(NAXES(1)*NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL AVERAGE(FLATFILE,(/FPIXELS(1),FPIXELS(2),1/),
     &  (/LPIXELS(1),LPIXELS(2),NAXES(3)/),DFLAT)
      CALL FTOPEN(UNIT,DARKFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(DARKFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      IF ((NAXES(1).NE.M) .OR. (NAXES(2).NE.N))THEN
        PRINT *,TRIM(TARGFILE)//' and '//TRIM(DARKFILE)//
     &    ' differ in size.'
        RETURN
      END IF
      CALL FTCLOS(UNIT, STATUS)
      ALLOCATE(DDARK(NAXES(1)*NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL AVERAGE(DARKFILE,(/FPIXELS(1),FPIXELS(2),1/),
     &  (/LPIXELS(1),LPIXELS(2),NAXES(3)/),DDARK)
      IF (MAXVAL(DDARK) .GE. MINVAL(DFLAT))THEN
        PRINT *,'Maximum of dark current larger than or equal'
     &    //' to minimum of flat field.'
        RETURN
      END IF
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=M*N
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      ALLOCATE(DTARG(NPIXELS))
      ALLOCATE(DCRCT(NPIXELS))
      CALL DELETEFILE(CRCTFILE,STATUS)
      DO K=FPIXELS(3),LPIXELS(3)
        CALL READIMAGE(TARGFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &    (/LPIXELS(1),LPIXELS(2),K/),DTARG)
        DCRCT=(DTARG-DDARK)/(DFLAT-DDARK)
        CALL APPENDIMAGE(CRCTFILE,(/1,1,K/),(/M,N,K/),DCRCT)
      END DO
      DEALLOCATE(DTARG)
      DEALLOCATE(DFLAT)
      DEALLOCATE(DDARK)
      DEALLOCATE(DCRCT)
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE FLATFIELDCORRECT
C ******************************************************************************
      SUBROUTINE ESTSNR(IMGFILE,FLAGFILE,FITN,DSNR,PREFIX)
C  Estimate the signal-to-noise ratio of given image.
C
C  Purpose:
C  ========
C  estimate_snr = var(signal) / var(noise)
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FITN
      INTEGER :: STATUS,NAXES(3),K,X,Y,NSAMPLES
      INTEGER, PARAMETER :: NBIN=20
      DOUBLE PRECISION, INTENT(OUT) :: DSNR
      DOUBLE PRECISION :: DMAX,DMIN,DMU,DSIGMA,DSVAR,DNVAR
      DOUBLE PRECISION, ALLOCATABLE :: DIMG(:,:),DFLAG(:,:),DBG(:,:),
     &  DB(:),DHIST(:)
      CHARACTER(LEN=*), INTENT(IN) :: IMGFILE,FLAGFILE,PREFIX
      INTERFACE
      SUBROUTINE BGFIT2PN(NX,NY,DFLAG,DIMG,N,DBG,DB)
      INTEGER, INTENT(IN) :: NX,NY,N
      DOUBLE PRECISION, INTENT(IN) :: DFLAG(NX,NY),DIMG(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(NX,NY),DB(*)
      END SUBROUTINE BGFIT2PN
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(IN) :: NAXES(3)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE IMAGESIZE
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
      END INTERFACE
      CALL IMAGESIZE(IMGFILE,NAXES)
      WRITE(*,'(A,I3,A,I3)')' input image size: ',NAXES(1),' x ',
     &  NAXES(2)
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(IMGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      ALLOCATE(DFLAG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(FLAGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DFLAG)
      NSAMPLES=NINT(SUM(DFLAG))
      WRITE(*,'(A,I5)')' number of flagged samples:',NSAMPLES
      ALLOCATE(DBG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DB(NAXES(1)*NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL BGFIT2PN(NAXES(1),NAXES(2),DFLAG,DIMG,FITN,DBG,DB)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_snr_bg.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'fitted background: '//TRIM(PREFIX)//'_snr_bg.fits'
      ALLOCATE(DHIST(NBIN),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      DBG=(DBG-DIMG)*DFLAG
      CALL WRITEIMAGE(TRIM(PREFIX)//'_snr_res.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'background residual: '//TRIM(PREFIX)//'_snr_res.fits'
      DMAX=MAXVAL(DBG)
      DMIN=MINVAL(DBG)
      DHIST=0.0D0
      DO X=1,NAXES(1)
        DO Y=1,NAXES(2)
          DHIST(NINT((DBG(X,Y)-DMIN)/(DMAX-DMIN)*DBLE(NBIN-1))+1)=1.0D0+
     &      DHIST(NINT((DBG(X,Y)-DMIN)/(DMAX-DMIN)*DBLE(NBIN-1))+1)
        END DO
      END DO
      PRINT *,'parameters of histogram:'
      PRINT *,'========================'
      DHIST=DHIST/DBLE(NSAMPLES)
      DO K=1,NBIN
        WRITE(*,'(A,I2,A,ES10.3)')' x=',K,', y=',DHIST(K)
      END DO
      PRINT *,'========================'
      DMU=SUM(DBG)/DBLE(NSAMPLES)
      DBG=DBG-DMU
      DNVAR=SUM(DBG*DBG)/DBLE(NSAMPLES)
      DSIGMA=DSQRT(DNVAR)
      DBG=DBG/DSIGMA
      WRITE(*,'(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)')
     &  ' min = ',DMIN,', max = ',DMAX,', mean = ',DMU,', std = ',DSIGMA
      DMU=SUM(DIMG)/DBLE(NAXES(1)*NAXES(2))
      WRITE(*,'(A,ES9.2)')' mean value of signal (given image): ',DMU
      DIMG=DIMG-DMU
      DSVAR=SUM(DIMG*DIMG)/DBLE(NAXES(1)*NAXES(2))
      WRITE(*,'(A,ES9.2,ES9.2)')' variance of signal and noise: ',
     &  DSVAR,DNVAR
      DSNR=DSVAR/DNVAR
      WRITE(*,'(A,ES9.2)')' estimated SNR: ',DSNR
      DEALLOCATE(DIMG)
      DEALLOCATE(DFLAG)
      DEALLOCATE(DBG)
      DEALLOCATE(DB)
      DEALLOCATE(DHIST)
      RETURN
      END SUBROUTINE ESTSNR
C ******************************************************************************
      SUBROUTINE BGFIT2PN(NX,NY,DFLAG,DIMG,N,DBG,DB)
C  2-dimensional background fitting subroutine using N-th polynomials.
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX,NY,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS,X,Y,L,PX,PY,K,
     &  STATUS
      DOUBLE PRECISION, INTENT(IN) :: DIMG(NX,NY),DFLAG(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(NX,NY),DB(*)
      DOUBLE PRECISION :: DXC,DYC
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:),DXI(:,:),DYI(:,:),
     &  DSAMPLEX(:),DSAMPLEY(:)
      STATUS=0
      DXC=0.5D0*(DBLE(NX+1))
      DYC=0.5D0*(DBLE(NY+1))
      NPARAMS=(N+1)*(N+2)/2
      NSAMPLES=NINT(SUM(DFLAG))
      WRITE(*,'(I6,A)'),NSAMPLES,' pixels have been sampled to'//
     &  ' determine the fitted background.'
      ALLOCATE(DSAMPLEX(NSAMPLES),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DSAMPLEY(NSAMPLES),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DXI(NX,NY),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DYI(NX,NY),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      NSAMPLES=0
      DO X=1,NX
        DO Y=1,NY
          DXI(X,Y)=DBLE(X)
          DYI(X,Y)=DBLE(Y)
          IF(DFLAG(X,Y) .EQ. 1.0D0)THEN
            NSAMPLES=NSAMPLES+1
            DB(NSAMPLES)=DIMG(X,Y)
            DSAMPLEX(NSAMPLES)=DBLE(X)
            DSAMPLEY(NSAMPLES)=DBLE(Y)
          END IF
        END DO
      END DO
      IF(N.EQ.-1)THEN
        PRINT *,'use maximum value.'
        DB(1)=MAXVAL(DB(1:NSAMPLES))
        DBG=DB(1)
        RETURN
      END IF
      IF(N.EQ.0)THEN
        PRINT *,'use constant average value.'
      ELSE
        WRITE(*,'(A,I2,A)')' use ',N,'-th polynomials.'
      END IF
      ALLOCATE(A(NSAMPLES,NPARAMS),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      L=1
      DO K=0,N
        DO PX=0,K
          PY=K-PX
          A(1:NSAMPLES,L)=DEXP(DLOG(DSAMPLEX(1:NSAMPLES))*DBLE(PX))*
     &      DEXP(DLOG(DSAMPLEY(1:NSAMPLES))*DBLE(PY))
          L=L+1
        END DO
      END DO
      DEALLOCATE(DSAMPLEX)
      DEALLOCATE(DSAMPLEY)
      LDA=NSAMPLES
      LDB=MAX(NSAMPLES,NPARAMS)
      LWORK=-1
      ALLOCATE(WORK(1),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,DB,LDB,WORK,
     &  LWORK,INFO)
      LWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,DB,LDB,WORK,
     &  LWORK,INFO)
      DEALLOCATE(WORK)
      DEALLOCATE(A)
      IF (INFO.GT.0)THEN
        PRINT *,'The least-squares solution could not be computed ',
     &    'because A does not have full rank.'
        RETURN
      ENDIF
      IF (INFO.LT.0)THEN
        PRINT *,'The argument has illegal value: ',ABS(INFO)
        RETURN
      ENDIF
      DBG=0.0D0
      L=1
      DO K=0,N
        DO PX=0,K
          PY=K-PX
          DBG=DBG+DB(L)*DEXP(DLOG(DXI)*DBLE(PX))*
     &      DEXP(DLOG(DYI)*DBLE(PY))
          L=L+1
        END DO
      END DO
      DEALLOCATE(DXI)
      DEALLOCATE(DYI)
      RETURN
      END SUBROUTINE BGFIT2PN
C ******************************************************************************
      SUBROUTINE BGFIT2P0(M,N,D,IMG,BG,B)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      NPARAMS=1
      NSAMPLES=0
      B(1)=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
            B(1)=B(1)+IMG(I,J)
          END IF
        END DO
      END DO
      B(1)=B(1)/DBLE(NSAMPLES)
      BG=B(1)
      RETURN
      END SUBROUTINE BGFIT2P0
C ******************************************************************************
      SUBROUTINE BGFIT2P2(M,N,D,IMG,BG,B)
C  2-dimensional background fitting subroutine using 2nd polynomials.
C
C  Purpose
C  =======
C  z = a_0 + a_1*x + a_2*y + a_3*x^2 + a_4*x*y + a_5*y^2
C  x is column number.
C  y is row number.
C  z is value of IMG at (x,y).
C  Try to determine parameters a_i, i=0,1,2,3,4,5.
C
C  Arguments
C  =========
C  M is the number of rows of matrix IMG.
C  N is the number of columns of matrix IMG.
C  IMG is the matrix the background of which is to be fitted.
C  BG,B are the output of this subroutine.
C  Each row of A is vector (1, x, y, x^2, x*y, y^2) for a specific (x,y).
C  B is (z_1, z_2, z_3, ..., z_n)'. The subscript of z denotes different
C  location.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:)
      NPARAMS=6
      NSAMPLES=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(A(NSAMPLES,6))
      K=1
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            A(K,:)=(/DBLE(1),X,Y,X*X,X*Y,Y*Y/)
            B(K)=IMG(INT(Y),INT(X))
            K=K+1
          ENDIF
        ENDDO
      ENDDO
      LDA=NSAMPLES
      LDB=MAX(NSAMPLES,6)
      LWORK=-1
      ALLOCATE(WORK(1))
      CALL DGELS('N',NSAMPLES,6,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      LWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DGELS('N',NSAMPLES,6,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      IF (INFO.GT.0)THEN
        PRINT *,'The least-squares solution could not be computed ',
     &    'because A does not have full rank.'
        RETURN
      ENDIF
      IF (INFO.LT.0)THEN
        PRINT *,'The argument has illegal value: ',ABS(INFO)
        RETURN
      ENDIF
      DO I=1,M
        DO J=1,N
          BG(I,J)=B(1)+B(2)*DBLE(J)+B(3)*DBLE(I)+B(4)*DBLE(J)*DBLE(J)
     &      +B(5)*DBLE(J)*DBLE(I)+B(6)*DBLE(I)*DBLE(I)
        ENDDO
      ENDDO
      DEALLOCATE(WORK)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE BGFIT2P2
C ******************************************************************************
      SUBROUTINE BGFIT2P4(M,N,D,IMG,BG,B)
C  2-dimensional background fitting subroutine using 4th polynomials.
C
C  Purpose
C  =======
C  z = a_0 + 
C    a_1*x + a_2*y + 
C    a_3*x^2 + a_4*x*y + a_5*y^2 + 

C    a_6*x^3 + a_7*x^2*y + a_8*x*y^2 + a_9*y^3 +
C    a_10*x^4 + a_11*x^3*y + a_12*x^2*y^2 + a_13*x*y^3 + a_14*y^4
C  x is column number.
C  y is row number.
C  z is value of IMG at (x,y).
C  Try to determine parameters a_i, i=0,1,2,3,...,14.
C
C  Arguments
C  =========
C  M is the number of rows of matrix IMG.
C  N is the number of columns of matrix IMG.
C  D is distance. Elements whose distance from the centre of IMG is larger than
C    D are sampled. BG is fitted to the sampled elements.
C  IMG is the matrix the background of which is to be fitted.
C  BG is the output of this subroutine.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:)
      NPARAMS=15
      NSAMPLES=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(A(NSAMPLES,NPARAMS))
      K=1
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (DSQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            A(K,:)=(/DBLE(1),X,Y,X*X,X*Y,Y*Y,
     &        X*X*X,X*X*Y,X*Y*Y,Y*Y*Y,
     &        X*X*X*X,X*X*X*Y,X*X*Y*Y,X*Y*Y*Y,Y*Y*Y*Y/)
            B(K)=IMG(INT(Y),INT(X))
            K=K+1
          ENDIF
        ENDDO
      ENDDO
      LDA=NSAMPLES
      LDB=MAX(NSAMPLES,NPARAMS)
      LWORK=-1
      ALLOCATE(WORK(1))
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      LWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      IF (INFO.GT.0)THEN
        PRINT *,'The least-squares solution could not be computed ',
     &    'because A does not have full rank.'
        RETURN
      ENDIF
      IF (INFO.LT.0)THEN
        PRINT *,'The argument has illegal value: ',ABS(INFO)
        RETURN
      ENDIF
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          BG(I,J)=B(1)+B(2)*X+B(3)*Y+
     &      B(4)*X*X+B(5)*X*Y+B(6)*Y*Y+
     &      B(7)*X*X*X+B(8)*X*X*Y+B(9)*X*Y*Y+B(10)*Y*Y*Y+
     &      B(11)*X*X*X*X+B(12)*X*X*X*Y+B(13)*X*X*Y*Y+B(14)*X*Y*Y*Y+
     &      B(15)*Y*Y*Y*Y
        ENDDO
      ENDDO
      DEALLOCATE(WORK)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE BGFIT2P4
C ******************************************************************************
      SUBROUTINE GETSNR(M,N,D,SNR,IMG,FTMETHOD)
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN1D,PLAN2D
      INTEGER :: M,N,I,J,K,NPIXS,NSPLS,INFO
      DOUBLE PRECISION :: SNR,D,X,Y,XC,YC,PS,PN
      DOUBLE PRECISION :: IMG(M,N),CLN(M,N),WORK(M,N)
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:)
      DOUBLE COMPLEX :: ZIN2D(M,N),ZOUT2D(M,N)
      DOUBLE COMPLEX, ALLOCATABLE :: ZIN1D(:),ZOUT1D(:)
      CHARACTER*(*) :: FTMETHOD
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
        RETURN
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      NPIXS=M*N
      NSPLS=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSPLS=NSPLS+1
          END IF
        END DO
      END DO
      ALLOCATE(BUFFER(NSPLS))
      ALLOCATE(ZIN1D(NSPLS))
      ALLOCATE(ZOUT1D(NSPLS))
      CALL DFFTW_PLAN_DFT_1D(PLAN1D,NSPLS,ZIN1D,ZOUT1D,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLAN2D,M,N,ZIN2D,ZOUT2D,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      IF (FTMETHOD .EQ. 'P2') THEN
        CALL BGFIT2P2(M,N,D,IMG,WORK,BUFFER)
      ELSE IF (FTMETHOD .EQ. 'P4') THEN
        CALL BGFIT2P4(M,N,D,IMG,WORK,BUFFER)
      ELSE IF (FTMETHOD .EQ. 'P0') THEN
        CALL BGFIT2P0(M,N,D,IMG,WORK,BUFFER(1))
      ELSE
        PRINT *,' Unknown fitting method.'
        RETURN
      END IF
      CLN=IMG-WORK
      K=0
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            K=K+1
            BUFFER(K)=CLN(I,J)
          END IF
        END DO
      END DO
      ZIN2D(1:M,1:N)=DCMPLX(CLN(1:M,1:N))
      CALL DFFTW_EXECUTE_DFT(PLAN2D,ZIN2D,ZOUT2D)
      PS=SUM(ZOUT2D*CONJG(ZOUT2D))/DBLE(NPIXS)/DSQRT(DBLE(NPIXS))
      ZIN1D(1:NSPLS)=DCMPLX(BUFFER(1:NSPLS))
      CALL DFFTW_EXECUTE_DFT(PLAN1D,ZIN1D,ZOUT1D)
      PN=SUM(ZOUT1D*CONJG(ZOUT1D))/DBLE(NSPLS)/DSQRT(DBLE(NSPLS))
      SNR=PS/PN
      CALL DFFTW_DESTROY_PLAN(PLAN1D)
      CALL DFFTW_DESTROY_PLAN(PLAN2D)
C     CALL DFFTW_CLEANUP_THREADS()
      DEALLOCATE(ZIN1D)
      DEALLOCATE(BUFFER)
      DEALLOCATE(ZOUT1D)
      RETURN
      END SUBROUTINE GETSNR
C ******************************************************************************

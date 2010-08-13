      SUBROUTINE ISAREC(OBSFILE,FPIXELS,LPIXELS,M,N,DREF,DBG,DISA,DR,
     &  NUMIT,PREFIX)
C
C  Arguments:
C  ==========
C  OBSFILE - Observed fits file to process.
C  FPIXELS - First pixels in each frame of the observed fits file.
C  LPIXELS - Last pixels in each frame of the observed fits file.
C  M       - Number of rows of DREF, DBG, and DISA.
C  N       - Number of columns of DREF, DBG, and DISA.
C  DREF    - Image of the reference star.
C  DBG     - Background of the average of all frames of the observed fits file.
C  DISA    - Input/Output. As input it is the initial estimate.
C  DR      - Radius of the border of the signal in the image.
C  NUMIT   - Number of iterations.
C  PREFIX  - Prefix of output filename.
C
C
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLANF,PLANB
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3),M,N,NUMIT
      INTEGER :: NPIXS,NUMFRM,K,L,INFO,X,Y
      DOUBLE PRECISION, INTENT(IN) :: DREF(M,N),DBG(M,N),DR
      DOUBLE PRECISION, INTENT(INOUT) :: DISA(M,N)
      DOUBLE PRECISION :: DSNR,DRMS,DBETA
      DOUBLE PRECISION, DIMENSION(M,N) :: DIMG,DCORE,WORK,DCORR
      DOUBLE COMPLEX,DIMENSION(M,N) :: ZISA,ZIN,ZOUT
      CHARACTER*(*), INTENT(IN) :: OBSFILE,PREFIX
      CHARACTER(LEN=256) :: ISAFILE,NUMSTR
      NPIXS=M*N
      NUMFRM=LPIXELS(3)-FPIXELS(3)+1
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_PLAN_DFT_2D(PLANF,M,N,ZIN,ZOUT,-1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLANB,M,N,ZIN,ZOUT,1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      DO L=1,NUMIT
        WRITE(*,'(A,I4)') 'Loop =',L
        WRITE(NUMSTR,*) L
C       CALL GETSNR(M,N,DR,DSNR,DISA,0)
        DBETA=0.5
        CALL DECONVCLEAN(M,N,DISA,DCORE,DREF,DBETA,100000)
        ISAFILE=TRIM(PREFIX)//'_ISA_CORE_'//
     &    TRIM(ADJUSTL(NUMSTR))//'.FITS'
        CALL WRITEIMAGE(ISAFILE,(/1,1,1/),(/M,N,1/),DCORE)
        WORK=0.0
        ZIN=CMPLX(DCORE)
        CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
        ZISA=CONJG(ZOUT)
        DO K=FPIXELS(3),LPIXELS(3)
          CALL READIMAGE(OBSFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &      (/LPIXELS(1),LPIXELS(2),K/),DIMG)
          DIMG=DIMG/SUM(DIMG)*DBLE(NPIXS)-DBG
          DIMG=DIMG/SUM(DIMG)*DBLE(NPIXS)
          ZIN=CMPLX(DIMG)
          CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
          ZIN=ZOUT*ZISA
          CALL DFFTW_EXECUTE_DFT(PLANB,ZIN,ZOUT)
          DCORR=DBLE(ZOUT)
C         CALL WRITEIMAGE('CORR.FITS',(/1,1,1/),(/M,N,1/),DCORR)
          X=MAXLOC(MAXVAL(DCORR,1),2)
          Y=MAXLOC(MAXVAL(DCORR,2),1)
          IF (N-X.GT.X-1)THEN
            X=X-1
          ELSE
            X=X-N
          END IF
          IF (M-Y.GT.Y-1)THEN
            Y=Y-1
          ELSE
            Y=Y-M
          END IF
          WORK=WORK+EOSHIFT(EOSHIFT(DIMG,Y-1,DBLE(0),1),
     &      X-1,DBLE(0),2)
        END DO
        WORK=WORK/DBLE(NUMFRM)
        ISAFILE=TRIM(PREFIX)//'_ISA_'//TRIM(ADJUSTL(NUMSTR))//'.FITS'
        CALL WRITEIMAGE(ISAFILE,(/1,1,1/),(/M,N,1/),WORK)
        DRMS=SQRT(SUM((WORK-DISA)*(WORK-DISA))/DBLE(NPIXS))
        WRITE(*,'(A,A,ES10.3)') TRIM(ISAFILE),', RMS =',DRMS
        DISA=WORK
      END DO
      CALL DFFTW_DESTROY_PLAN(PLANF)
      CALL DFFTW_DESTROY_PLAN(PLANB)
      RETURN
      END SUBROUTINE ISAREC

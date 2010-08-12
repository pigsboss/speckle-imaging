      SUBROUTINE ISAREC(OBSFILE,FPIXELS,LPIXELS,M,N,DREF,DBG,DISA,DR,
     &  P,NUMIT,PREFIX)
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
C  P       - Pad all matrices to P x P for fast fourier transform.
C  NUMIT   - Number of iterations.
C  PREFIX  - Prefix of output filename.
C
C
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN,PLANI
      INTEGER :: FPIXELS(*),LPIXELS(*),NUMIT,K,L,INFO,X,Y
      INTEGER :: M,N,NPIXS,NUMFRM,P
      DOUBLE PRECISION :: DSNR,DR,DRMS
      DOUBLE PRECISION, DIMENSION(M,N) :: DIMG,DREF,DBG,DISA,DCORE
      DOUBLE PRECISION, DIMENSION(M,N) :: WORK
      DOUBLE PRECISION, DIMENSION(P,P) :: DCORR
      DOUBLE COMPLEX,DIMENSION(P,P) :: ZISA,ZIN,ZOUT
      CHARACTER*(*) :: PREFIX
      CHARACTER*(256) :: OBSFILE,ISAFILE,NUMSTR
      NPIXS=M*N
      NUMFRM=LPIXELS(3)-FPIXELS(3)+1
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_PLAN_DFT_2D(PLAN,P,P,ZIN,ZOUT,-1,
     &  FFTW_PATIENT+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLANI,P,P,ZIN,ZOUT,1,
     &  FFTW_PATIENT+FFTW_DESTROY_INPUT)
      DO L=1,NUMIT
        WRITE(*,'(A,I4)') 'Loop =',L
        WRITE(NUMSTR,*) L
        CALL GETSNR(M,N,DR,DSNR,DISA,0)
        CALL DECONVWNR(M,N,DISA,DCORE,DREF,DSNR)
        ISAFILE=TRIM(PREFIX)//'_ISA_CORE_'//
     &    TRIM(ADJUSTL(NUMSTR))//'.FITS'
        CALL WRITEIMAGE(ISAFILE,(/1,1,1/),(/M,N,1/),DCORE)
        WORK=0.0
        ZIN=CMPLX(0.0)
        ZIN(1:M,1:N)=CMPLX(DCORE)
        CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
        ZISA=CONJG(ZOUT)
        DO K=FPIXELS(3),LPIXELS(3)
          CALL READIMAGE(OBSFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &      (/LPIXELS(1),LPIXELS(2),K/),DIMG)
          DIMG=DIMG/SUM(DIMG)*DBLE(NPIXS)-DBG
          DIMG=DIMG/SUM(DIMG)*DBLE(NPIXS)
          ZIN=CMPLX(0.0)
          ZIN(1:M,1:N)=CMPLX(DIMG)
          CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
          ZIN=ZOUT*ZISA
          CALL DFFTW_EXECUTE_DFT(PLANI,ZIN,ZOUT)
          DCORR=DBLE(ZOUT)
C         CALL WRITEIMAGE('CORR.FITS',(/1,1,1/),(/P,P,1/),DCORR)
          X=MAXLOC(MAXVAL(DCORR,1),2)
          Y=MAXLOC(MAXVAL(DCORR,2),1)
          IF (P-X.GT.X-1)THEN
            X=X-1
          ELSE
            X=X-P
          END IF
          IF (P-Y.GT.Y-1)THEN
            Y=Y-1
          ELSE
            Y=Y-P
          END IF
          WORK=WORK+EOSHIFT(EOSHIFT(DIMG,Y-1,DBLE(0),1),
     &      X-1,DBLE(0),2)
        END DO
        WORK=WORK/DBLE(NUMFRM)
        ISAFILE=TRIM(PREFIX)//'_ISA_'//TRIM(ADJUSTL(NUMSTR))//'.FITS'
        CALL WRITEIMAGE(ISAFILE,(/1,1,1/),(/M,N,1/),WORK)
        DRMS=SQRT(SUM((WORK-DISA)*(WORK-DISA))/DBLE(NPIXS))
        WRITE(*,'(A,A,ES10.3,A,ES10.3)') TRIM(ISAFILE),', RMS =',DRMS,
     &    ' SNR =',DSNR
        DISA=WORK
      END DO
      CALL DFFTW_DESTROY_PLAN(PLAN)
      CALL DFFTW_DESTROY_PLAN(PLANI)
      CALL DFFTW_CLEANUP_THREADS()
      RETURN
      END SUBROUTINE ISAREC
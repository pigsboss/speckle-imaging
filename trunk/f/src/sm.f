      PROGRAM MAIN
C  Speckle masking main program.
C
C  Usage:
C  ======
C  sm file_obs fpixels_obs lpixels_obs radius_obs
C    file_ref fpixels_ref lpixels_ref radius_ref
C    fit_method max_width prefix
C
C  Declaration:
C  ============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER :: FPOBS(3),LPOBS(3),FPREF(3),LPREF(3),
     &  P,MW,BUFFERSIZE,M,N,NBISP,K,L
      INTEGER*8 :: PLAN
      DOUBLE PRECISION :: DROBS,DRREF,DSNR,DTMP
      DOUBLE PRECISION, ALLOCATABLE :: DAVG(:,:),DBETA(:),DPHI(:,:)
     &  ,DRHO(:,:),DPSDOBS(:,:),DPSDREF(:,:),DSM(:,:)
      DOUBLE COMPLEX :: ZTMP
      DOUBLE COMPLEX, ALLOCATABLE :: ZBISPOBS(:),ZBISPREF(:),
     & ZIN(:,:),ZOUT(:,:)
      CHARACTER(LEN=256) :: FILEOBS,ARG,FILEREF,PREFIX
      CHARACTER(LEN=10) :: FTMETHOD
C
C  Statements:
C  ===========
C
      CALL GETARG(1,FILEOBS)
      CALL GETARG(2,ARG)
      READ(ARG,*) FPOBS(1),FPOBS(2),FPOBS(3)
      CALL GETARG(3,ARG)
      READ(ARG,*) LPOBS(1),LPOBS(2),LPOBS(3)
      CALL GETARG(4,ARG)
      READ(ARG,*) DROBS
      CALL GETARG(5,FILEREF)
      CALL GETARG(6,ARG)
      READ(ARG,*) FPREF(1),FPREF(2),FPREF(3)
      CALL GETARG(7,ARG)
      READ(ARG,*) LPREF(1),LPREF(2),LPREF(3)
      CALL GETARG(8,ARG)
      READ(ARG,*) DRREF
      CALL GETARG(9,FTMETHOD)
      CALL GETARG(10,ARG)
      READ(ARG,*) MW
      CALL GETARG(11,PREFIX)
C  The padding size must be even:
      M=LPOBS(1)-FPOBS(1)+1
      N=LPOBS(2)-FPOBS(2)+1
      P=MAX(M,N)
      DTMP=DEXP(CEILING(DLOG(DBLE(P))/DLOG(2.D0))*DLOG(2.D0))
      IF (DTMP-FLOOR(DTMP) .GE. 0.D5)THEN
        P=INT(FLOOR(DTMP)+1)
      ELSE
        P=INT(FLOOR(DTMP))
      END IF
      NBISP=(MW+1)*(2*P-MW)*P*(P+2)/8
      PRINT *,'Padding size: ',P
      PRINT *,'Bispectrum size (MB): ',DBLE(NBISP)*16/1024/1024
      ALLOCATE(DAVG(M,N))
      CALL AVERAGE(FILEOBS,FPOBS,LPOBS,DAVG,10)
      CALL GETSNR(M,N,DROBS,DSNR,DAVG,FTMETHOD)
      ALLOCATE(ZBISPOBS(NBISP))
      CALL BISPECTRUM(FILEOBS,FPOBS,LPOBS,DAVG,DROBS,FTMETHOD,
     &  P,MW,ZBISPOBS)
      M=LPREF(1)-FPREF(1)+1
      N=LPREF(2)-FPREF(2)+1
      DEALLOCATE(DAVG)
      ALLOCATE(DAVG(M,N))
      CALL AVERAGE(FILEREF,FPREF,LPREF,DAVG,10)
      ALLOCATE(ZBISPREF(NBISP))
      CALL BISPECTRUM(FILEREF,FPREF,LPREF,DAVG,DRREF,FTMETHOD,
     &  P,MW,ZBISPREF)
      ALLOCATE(DBETA(NBISP))
      DO K=1,NBISP
        IF (ZABS(ZBISPREF(K))>DBLE(1.0)/DSNR)THEN
          ZTMP=ZBISPOBS(K)/ZBISPREF(K)
          DBETA(K)=DATAN2(DIMAG(ZTMP),DREAL(ZTMP))
        ELSE
          DBETA(K)=0
        END IF
      END DO
      DEALLOCATE(ZBISPOBS)
      DEALLOCATE(ZBISPREF)
      ALLOCATE(DPHI(P,P))
      CALL PHASERECURSION(P,MW,DBETA,DPHI)
      DEALLOCATE(DBETA)
      ALLOCATE(DPSDOBS(P,P))
      ALLOCATE(DPSDREF(P,P))
      CALL GETPSD(FILEOBS,FPOBS,LPOBS,DROBS,P,FTMETHOD,DPSDOBS)
      CALL GETPSD(FILEREF,FPREF,LPREF,DRREF,P,FTMETHOD,DPSDREF)
      ALLOCATE(DRHO(P,P))
      DRHO=DPSDOBS/DPSDREF
      ALLOCATE(ZIN(P,P))
      ALLOCATE(ZOUT(P,P))
      CALL DFFTW_PLAN_DFT_2D(PLAN,P,P,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=CMPLX(DRHO*DCOS(DPHI),DRHO*DSIN(DPHI))
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_RHO.FITS',(/1,1,1/),(/P,P,1/),
     &  DRHO)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_PHI.FITS',(/1,1,1/),(/P,P,1/),
     &  DPHI)
      ALLOCATE(DSM(P,P))
      DSM=DREAL(ZOUT)
      CALL DFFTSHIFT(P,P,DSM)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_SMK.FITS',(/1,1,1/),(/P,P,1/),
     &  DSM)
      DEALLOCATE(DSM)
      DEALLOCATE(DPHI)
      DEALLOCATE(ZIN)
      DEALLOCATE(ZOUT)
      DEALLOCATE(DRHO)
      DEALLOCATE(DPSDOBS)
      DEALLOCATE(DPSDREF)
      CALL DFFTW_DESTROY_PLAN(PLAN)
      STOP
      END PROGRAM MAIN
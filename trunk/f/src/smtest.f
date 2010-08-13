      PROGRAM SMTEST
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER :: N,NFRAMES,MW,NBISP
      INTEGER*8 :: PLAN
      DOUBLE PRECISION, ALLOCATABLE :: DBETA(:),DPHI(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZBISP(:),ZIN(:,:),ZOUT(:,:)
      CHARACTER(LEN=256) :: ARG,FTMETHOD
      CALL GETARG(1,ARG)
      READ(ARG,*) N
      CALL GETARG(2,ARG)
      READ(ARG,*) NFRAMES
      CALL GETARG(3,ARG)
      READ(ARG,*) MW
      NBISP=(MW+1)*(2*N-MW)*N*(N+2)/8
      FTMETHOD='P0'
      ALLOCATE(ZBISP(NBISP))
      CALL BISPECTRUM('SMTEST.FITS',(/1,1,1/),(/N,N,NFRAMES/),MW,ZBISP)
      ALLOCATE(DBETA(NBISP))
      DBETA=DATAN2(DIMAG(ZBISP),DREAL(ZBISP))
      DEALLOCATE(ZBISP)
      ALLOCATE(DPHI(N,N))
      CALL PHASERECURSION(N,N,MW,DBETA,DPHI)
      ALLOCATE(ZIN(N,N))
      ALLOCATE(ZOUT(N,N))
      CALL DFFTW_PLAN_DFT_2D(PLAN,N,N,ZIN,ZOUT,1,FFTW_ESTIMATE)
      CALL WRITEIMAGE('SMTESTPHASE.FITS',(/1,1,1/),(/N,N,1/),DPHI)
      ZIN=CMPLX(DCOS(DPHI),DSIN(DPHI))
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DEALLOCATE(DPHI)
      DEALLOCATE(DBETA)
      DEALLOCATE(ZIN)
      CALL ZFFTSHIFT(N,N,ZOUT)
      ZOUT=ZOUT/DBLE(N*N)
      CALL WRITEIMAGE('SMTESTOBJ.FITS',(/1,1,1/),(/N,N,1/),DREAL(ZOUT))
      DEALLOCATE(ZOUT)
      STOP
      END PROGRAM SMTEST
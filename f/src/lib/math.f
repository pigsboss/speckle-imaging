      SUBROUTINE POISSRNDM(NX,NY,DLAMBDA,RND)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NX,NY
      DOUBLE PRECISION,INTENT(IN) :: DLAMBDA(NX,NY)
      INTEGER,INTENT(OUT) :: RND(NX,NY)
      INTEGER :: K(NX,NY),FLAG(NX,NY)
      DOUBLE PRECISION :: DL(NX,NY),DP(NX,NY),DU(NX,NY)
      DL=DEXP(-1.0D0 * DLAMBDA)
      K=0
      DP=1.0D0
      FLAG=1
      DO
        K = K + FLAG
        CALL RANDOM_NUMBER(DU)
        DP = DP * DU
        FLAG = IDNINT(-1.0D0*DSIGN(0.5D0,DL-DP)+0.5D0)
        IF (SUM(FLAG) .LE. 0) THEN
          EXIT
        END IF
      END DO
      RND = K - 1
      RETURN
      END SUBROUTINE POISSRNDM
C ******************************************************************************
      FUNCTION POISSRND(DLAMBDA)
      IMPLICIT NONE
      INTEGER :: POISSRND
      DOUBLE PRECISION,INTENT(IN) :: DLAMBDA
      INTEGER :: K
      DOUBLE PRECISION :: DL,DP,DU
      DL=DEXP(-1.0D0 * DLAMBDA)
      K=0
      DP=1.0D0
      DO
        K = K + 1
        CALL RANDOM_NUMBER(DU)
        DP = DP * DU
        IF (DP .LE. DL) THEN
          EXIT
        END IF
      END DO
      POISSRND = K - 1
      RETURN
      END FUNCTION POISSRND
C ******************************************************************************
      SUBROUTINE INIT_RANDOM_SEED()
      IMPLICIT NONE
      INTEGER :: I, N, CLOCK
      INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
      CALL RANDOM_SEED(SIZE = N)
      ALLOCATE(SEED(N))
      CALL SYSTEM_CLOCK(COUNT = CLOCK)
      SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
      CALL RANDOM_SEED(PUT = SEED)
      DEALLOCATE(SEED)
      RETURN
      END SUBROUTINE INIT_RANDOM_SEED
C ******************************************************************************
      SUBROUTINE GETARGUMENT(NX,NY,ZX,DPHI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(OUT) :: DPHI(NX,NY)
      DOUBLE COMPLEX, INTENT(IN) :: ZX(NX,NY)
      INTEGER :: X,Y
      DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979323846D0
      DOUBLE PRECISION :: DX(NX,NY),DY(NX,NY),DTMP
      DX=DREAL(ZX)
      DY=DIMAG(ZX)
      DO X=1,NX
        DO Y=1,NY
          IF(DX(X,Y) .LT. 0.0D0 .AND. DY(X,Y) .LT. 0.0D0)THEN
            DTMP=DY(X,Y)/DX(X,Y)
            DPHI(X,Y)=-1.0D0*PI-DATAN(DTMP)
          ELSE IF(DX(X,Y) .EQ. 0.0D0 .AND. DY(X,Y) .LT. 0.0D0)THEN
            DPHI(X,Y)=-0.5D0*PI
          ELSE IF(DX(X,Y) .GT. 0.0D0 .AND. DY(X,Y) .LT. 0.0D0)THEN
            DTMP=DABS(DY(X,Y))/DX(X,Y)
            DPHI(X,Y)=0.0D0-DATAN(DTMP)
          ELSE IF(DX(X,Y) .GE. 0.0D0 .AND. DY(X,Y) .EQ. 0.0D0)THEN
            DPHI(X,Y)=0.0D0
          ELSE IF(DX(X,Y) .GT. 0.0D0 .AND. DY(X,Y) .GT. 0.0D0)THEN
            DTMP=DY(X,Y)/DX(X,Y)
            DPHI(X,Y)=DATAN(DTMP)
          ELSE IF(DX(X,Y) .EQ. 0.0D0 .AND. DY(X,Y) .GT. 0.0D0)THEN
            DPHI(X,Y)=0.5D0*PI
          ELSE IF(DX(X,Y) .LT. 0.0D0 .AND. DY(X,Y) .GT. 0.0D0)THEN
            DTMP=DY(X,Y)/DABS(DX(X,Y))
            DPHI(X,Y)=PI-DATAN(DTMP)
          ELSE IF(DX(X,Y) .LT. 0.0D0 .AND. DY(X,Y) .EQ. 0.0D0)THEN
            DPHI(X,Y)=PI
          ELSE
            WRITE(*,*)DX(X,Y),DY(X,Y)
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE GETARGUMENT
C ******************************************************************************
      SUBROUTINE MULTIPLYDP(NA,DPA,NB,DPB,DPC)
C  Multiplication of double precision polynomials.
C
C  Purpose:
C  ========
C  c = a x b
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NA,NB
      DOUBLE PRECISION, INTENT(IN) :: DPA(NA+1),DPB(NB+1)
      DOUBLE PRECISION, INTENT(OUT) :: DPC(NA+NB+1)
      INTEGER :: NC,KA,KB,KC
      NC=NB+NA
      DPC=0.0D0
      DO KA=0,NA
        DO KB=0,NB
          KC=KA+KB
          DPC(KC+1)=DPC(KC+1)+DPA(KA+1)*DPB(KB+1)
        END DO
      END DO
      RETURN
      END SUBROUTINE MULTIPLYDP
C ******************************************************************************
      SUBROUTINE INTEGRATEDP(NA,DPA,DXMIN,DXMAX,DS)
C  Definite integration of double precision polynomial function.
C
C  Purpose:
C  ========
C  s = \int_xmin^xmax a(x) dx
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NA
      INTEGER :: K
      DOUBLE PRECISION, INTENT(IN) :: DPA(NA+1),DXMIN,DXMAX
      DOUBLE PRECISION, INTENT(OUT) :: DS
      DOUBLE PRECISION :: DTMP
      DS=0.0D0
      DO K=0,NA
        DTMP=DEXP(DLOG(DXMAX)*DBLE(K))-DEXP(DLOG(DXMIN)*DBLE(K))
        DS=DS+DTMP*DPA(K+1)/DBLE(K+1)
      END DO
      RETURN
      END SUBROUTINE INTEGRATEDP
C ******************************************************************************
C      SUBROUTINE DLEGENDRE(DXMIN,DXMAX,NA,DLPA)
C      RETURN
C      END SUBROUTINE DLEGENDRE
C ******************************************************************************
      SUBROUTINE DCORR2D(M,N,DA,DB,DC)
C  correlation coefficients of 2-dimensional array.
C
C  Purpose:
C  ========
C  c = corr(a, b). calculates c given a and b using fast fourier transform.
C  fft(c) = fft(a) * conj(fft(b)) (periodic extension is implied).
C
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER :: M,N
      DOUBLE PRECISION, DIMENSION(M,N) :: DA,DB,DC
      DOUBLE COMPLEX, DIMENSION(M,N) :: ZA,ZB,ZIN,ZOUT
      CALL DFFTW_INIT_THREADS()
      CALL DFFTW_PLAN_WITH_NTHREADS(PLAN,2)
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=DCMPLX(DA)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZA=ZOUT
      ZIN=DCMPLX(DB)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZB=ZOUT
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZA*CONJG(ZB)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DC=DBLE(ZOUT)
      RETURN
      END SUBROUTINE DCORR2D

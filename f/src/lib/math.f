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
      END SUBROUTINE MULTIPLYP
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
      SUBROUTINE DLEGENDRE(DXMIN,DXMAX,NA,DLPA)
      RETURN
      END SUBROUTINE DLEGENDRE
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
      ZIN=CMPLX(DA)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZA=ZOUT
      ZIN=CMPLX(DB)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZB=ZOUT
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZA*CONJG(ZB)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DC=DBLE(ZOUT)
      RETURN
      END SUBROUTINE DCORR2D

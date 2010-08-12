      PROGRAM TESTSM
      IMPLICIT NONE
      INTEGER :: P,MW,I,J
      DOUBLE PRECISION, ALLOCATABLE :: DPHI(:,:),DBETA(:),DREF(:,:)
      CHARACTER :: ARG*256
      PRINT *,ZABS(COMPLEX(1.D0,1.D0))
      CALL GETARG(1,ARG)
      READ(ARG,*) P
      PRINT *,P
      CALL GETARG(2,ARG)
      READ(ARG,*) MW
      PRINT *,MW
      ALLOCATE(DPHI(P,P))
      ALLOCATE(DREF(P,P))
      ALLOCATE(DBETA((MW+1)*(2*P-MW)*P*(P+2)/8))
      DO I=1,P
        DO J=1,P
          DPHI(I,J)=DBLE(I)/DBLE(P)+(DBLE(J)/DBLE(P))*(DBLE(J)/DBLE(P))
        END DO
      END DO
      DPHI(1,1)=0
      DPHI(2,1)=0
      DPHI(1,2)=0
      DREF=DPHI
      CALL WRITEIMAGE('ORIGINAL.FITS',(/1,1,1/),(/P,P,1/),DPHI)
      CALL GENERATEBETA(P,MW,DPHI,DBETA)
      CALL PHASERECURSION(P,MW,DBETA,DPHI,DREF)
      CALL WRITEIMAGE('RECURSION.FITS',(/1,1,1/),(/P,P,1/),DPHI)
      DEALLOCATE(DPHI)
      DEALLOCATE(DBETA)
      STOP
      END PROGRAM TESTSM
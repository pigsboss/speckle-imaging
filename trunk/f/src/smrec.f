      SUBROUTINE BISPECTRUM(FILENAME,FPIXELS,LPIXELS,DAVG,DR,FTMETHOD,
     &  P,MW,ZBISP)
C  Calculate the mean bispectrum of all
C  the frames.
C
C  Arguments:
C  ==========
C  FILENAME   - Input filename.
C  FPIXELS    - First pixels.
C  LPIXELS    - Last pixels.
C  DAVG       - Previously computed average frame.
C  DR         - Radius of the border of signals in each image.
C  FTMETHOD   - Fitting method.
C  P          - Padding size.
C  MW         - Maximum width of the gaps between two frequency components.
C  BUFFERSIZE - Size of reading buffer in MB.
C  ZBISP      - Mean bispectrum.
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3),MW
      INTEGER, INTENT(IN) :: P
      INTEGER :: I1,I2,I3,J1,J2,K,LF,LT,LR,L,M,N,NPIXELS,NFRAMES,
     &  LBUFFER,INFO,BUFFERSIZE
      INTEGER*8 :: PLAN
      DOUBLE PRECISION, INTENT(IN) :: DR,
     &  DAVG(LPIXELS(1)-FPIXELS(1)+1,LPIXELS(2)-FPIXELS(2)+1)
      DOUBLE PRECISION, INTENT(OUT) :: ZBISP((MW+1)*(2*P-MW)*P*(P+2)/8)
      DOUBLE PRECISION, DIMENSION(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(2)-FPIXELS(2)+1) :: DBG,DIMG,DB
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:,:,:)
      DOUBLE COMPLEX, DIMENSION(0:P-1,0:P-1) :: ZIN,ZOUT,ZSP
      CHARACTER, INTENT(IN) :: FILENAME*256,FTMETHOD*10
      INTERFACE
      FUNCTION BISPOS(I1,J1,I2,J2,P)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1,J1,I2,J2,P
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      END INTERFACE
C
C  Statements:
C  ===========
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=M*N
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      BUFFERSIZE=10
C  Determine the length of the buffer:
      LBUFFER=INT(FLOOR(DBLE(BUFFERSIZE)*1024*1024/DBLE(8*NPIXELS)))
C  Calculate the background with given method:
      IF (TRIM(FTMETHOD).EQ.'P0')THEN
        CALL BGFIT2P0(M,N,DR,DAVG,DBG,DB(1,1))
      ELSE IF (TRIM(FTMETHOD).EQ.'P2')THEN
        CALL BGFIT2P2(M,N,DR,DAVG,DBG,DB)
      ELSE IF (TRIM(FTMETHOD).EQ.'P4')THEN
        CALL BGFIT2P4(M,N,DR,DAVG,DBG,DB)
      ELSE
        PRINT *,'Unknown fitting method ',TRIM(FTMETHOD)
        RETURN
      END IF
C  Calculate bispectrum:
      ALLOCATE(BUFFER(M,N,LBUFFER))
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
      PRINT *,'Start planning.'
      CALL DFFTW_PLAN_DFT_2D(PLAN,P,P,ZIN,ZOUT,-1,
     &  FFTW_PATIENT+FFTW_DESTROY_INPUT)
      PRINT *,'Finished planning.'
C  Determine the size of the bispectrum array:
      ZBISP=CMPLX(0)
      DIMG=0
      LF=1
      DO L=1,INT(CEILING(DBLE(NFRAMES)/DBLE(LBUFFER)))
        LT=MIN(LF+LBUFFER-1,LPIXELS(3))
        LR=LT-LF+1
        CALL READIMAGE(FILENAME,(/FPIXELS(1),FPIXELS(2),LF/),
     &    (/LPIXELS(1),LPIXELS(2),LT/),BUFFER(1:M,1:N,1:LBUFFER))
        DO K=1,LR
          DIMG=BUFFER(1:M,1:N,K)/SUM(BUFFER(1:M,1:N,K))*DBLE(NPIXELS)-
     &      DBG
          DIMG=DIMG/SUM(DIMG)*DBLE(NPIXELS)
          ZIN(0:M-1,0:N-1)=CMPLX(DIMG)
          CALL ZIFFTSHIFT(P,P,ZIN)
          CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
          ZSP=ZOUT
          I3=1
          DO J2=0,MW
            DO I2=0,P-1
              DO J1=0,P-1-J2
                DO I1=0,MIN(I2,P-1-I2)
C                 IF (I3 .NE. BISPOS(I1,J1,I2,J2,P))THEN
C                   PRINT *,I3,'(',I1,J1,I2,J2,')',BISPOS(I1,J1,I2,J2,P)
C                   RETURN
C                 END IF
                  ZBISP(I3)=ZBISP(I3)+ZSP(I1,J1)*ZSP(I2,J2)*
     &              CONJG(ZSP(I1+I2,J1+J2))
                  I3=I3+1
                END DO
              END DO
            END DO
          END DO
        END DO
        LF=LT+1
      END DO
      ZBISP=ZBISP/DBLE(NFRAMES)
      DEALLOCATE(BUFFER)
      RETURN
      END SUBROUTINE BISPECTRUM
C ******************************************************************************
      SUBROUTINE GENERATEBETA(P,MW,DPHI,DBETA)
C  Generate principal values of the arguments of the bispectrum of given phase
C  (beta) for testing the recursion algorithm.
C
      INTEGER, INTENT(IN) :: P,MW
      INTEGER :: I1,I2,J1,J2,K
      DOUBLE PRECISION, INTENT(IN) :: DPHI(0:P-1,0:P-1)
      DOUBLE PRECISION, INTENT(OUT) :: DBETA((MW+1)*(2*P-MW)*P*(P+2)/8)
C  Interface:
C  ==========
      INTERFACE
      FUNCTION BISPOS(I1,J1,I2,J2,P)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1,J1,I2,J2,P
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      END INTERFACE
      PRINT *,'GENERATE BETA'
      K=1
      DO J2=0,MW
        DO I2=0,P-1
          DO J1=0,P-1-J2
            DO I1=0,MIN(I2,P-1-I2)
              DBETA(K)=DPHI(I1,J1)+DPHI(I2,J2)-DPHI(I1+I2,J1+J2)
C             PRINT *,I1,J1,I2,J2,DPHI(I1,J1),DPHI(I2,J2),
C    &          DPHI(I1+I2,J1+J2),DBETA(K)
              K=K+1
            END DO
          END DO
        END DO
      END DO
      END SUBROUTINE GENERATEBETA
C ******************************************************************************
      SUBROUTINE PHASERECURSION(P,MW,DBETA,DPHI,DREF)
C  Recursive algorithm to solve the phase equations.
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: P,MW
      INTEGER :: I,J,K,I1,I2,J1,J2
      DOUBLE PRECISION :: R
      DOUBLE PRECISION, INTENT(OUT) :: DPHI(0:P-1,0:P-1)
      DOUBLE PRECISION, INTENT(IN) :: DBETA((MW+1)*(2*P-MW)*P*(P+2)/8)
     &  ,DREF(0:P-1,0:P-1)
C  Interface:
C  ==========
      INTERFACE
      FUNCTION BISPOS(I1,J1,I2,J2,P)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1,J1,I2,J2,P
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      END INTERFACE
C  Statements:
C  ===========
      PRINT *,'PHASE RECURSION'
      DPHI=DBLE(0.0)
      DO K=2,2*(P-1)
        DO I=MAX(0,K+1-P),MIN(K,P-1)
          J=K-I
          R=0.0
          DO I1=0,INT(FLOOR(DBLE(I)/2.0))
            I2=I-I1
            DO J2=0,MIN(J,MW)
              J1=J-J2
              IF (((I1.NE.I).OR.(J1.NE.J)).AND.((I2.NE.I).OR.(J2.NE.J)))
     &          THEN
                R=R+1.0
                DPHI(I,J)=DPHI(I,J)*(R-1.0)/R+
     &            (DPHI(I1,J1)+DPHI(I2,J2)-
     &            DBETA(BISPOS(I1,J1,I2,J2,P)))/R
              END IF
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE PHASERECURSION
C ******************************************************************************
      INTEGER FUNCTION BISPOS(I1,J1,I2,J2,P)
C  Calculate position in bispectrum array according to positions of the
C  frequency components in the spectrum matrices.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1,J1,I2,J2,P
      INTEGER :: K,L
      K=1+P*(P+2)*J2*(2*P-J2+1)/8
      IF (I2 .LT. P/2) THEN
        BISPOS=K+I2*(I2+1)*(P-J2)/2+J1*(I2+1)+I1
        RETURN
      END IF
      IF (I2 .EQ. P/2) THEN
        BISPOS=K+P*(P+2)*(P-J2)/8+J1*P/2+I1
        RETURN
      END IF
      IF (I2 .GT. P/2) THEN
        BISPOS=K+P*(P+2)*(P-J2)/8+
     &    (3*P-2*I2+2)*(2*I2-P)*(P-J2)/8+J1*(P-I2)+I1
        RETURN
      END IF
      END FUNCTION BISPOS
C *****************************************************************************

      SUBROUTINE BISPECTRUM(FILENAME,FPIXELS,LPIXELS,MW,DBG,ZBISP)
C  Calculate the mean bispectrum of all
C  the frames.
C
C  Arguments:
C  ==========
C  FILENAME   - Input filename.
C  FPIXELS    - First pixels.
C  LPIXELS    - Last pixels.
C  MW         - Maximum width of the gaps between two frequency components.
C  BUFFERSIZE - Size of reading buffer in MB.
C  ZBISP      - Mean bispectrum.
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3),MW
      INTEGER :: I1,I2,I3,J1,J2,K,LF,LT,LR,L,M,N,NPIXELS,NFRAMES,
     &  LBUFFER,BUFFERSIZE,LBISP,INFO
      INTEGER*8 :: PLAN
      DOUBLE COMPLEX, INTENT(OUT) :: ZBISP(*)
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:,:,:)
      DOUBLE PRECISION, INTENT(IN) :: DBG(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(1)-FPIXELS(1)+1)
      DOUBLE COMPLEX, DIMENSION(0:LPIXELS(1)-FPIXELS(1),
     &  0:LPIXELS(2)-FPIXELS(2)) :: ZIN,ZOUT,ZSP
      CHARACTER*(*), INTENT(IN) :: FILENAME
C
      INTERFACE
      FUNCTION BISPOS(I1,J1,I2,J2,M,N)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1,J1,I2,J2,M,N
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      END INTERFACE
C
C  Statements:
C  ===========
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      LBISP=(MW+1)*(2*N-MW)*M*(M+2)/8
      NPIXELS=M*N
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      BUFFERSIZE=8
C  Determine the length of the buffer:
      LBUFFER=INT(FLOOR(DBLE(BUFFERSIZE)*1024*1024/DBLE(8*NPIXELS)))
C  Calculate bispectrum:
      PRINT *,'Calculate the bispectrum.'
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
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,-1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      PRINT *,'Finished planning.'
      ZBISP(1:LBISP)=CMPLX(0)
      LF=1
      ALLOCATE(BUFFER(M,N,LBUFFER))
      DO L=1,INT(CEILING(DBLE(NFRAMES)/DBLE(LBUFFER)))
        LT=MIN(LF+LBUFFER-1,LPIXELS(3))
        LR=LT-LF+1
        CALL READIMAGE(FILENAME,(/FPIXELS(1),FPIXELS(2),LF/),
     &    (/LPIXELS(1),LPIXELS(2),LT/),BUFFER(1:M,1:N,1:LBUFFER))
        DO K=1,LR
          WRITE(*,'(A,I7,A,I7)') ' Processing ',K+LF-1,' of ',NFRAMES
          ZIN=CMPLX(BUFFER(1:M,1:N,K)-DBG)
          DO I1=0,M-1
            DO J1=0,N-1
              IF (ZABS(ZIN(I1,J1)).LT.DBLE(0))THEN
                ZIN(I1,J1)=CMPLX(0)
              END IF
            END DO
          END DO
          CALL ZIFFTSHIFT(M,N,ZIN)
          CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
          ZSP=ZOUT
          I3=1
          DO J2=0,MW
            DO I2=0,M-1
              DO J1=0,N-1-J2
                DO I1=0,MIN(I2,M-1-I2)
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
      DEALLOCATE(BUFFER)
      ZBISP(1:LBISP)=ZBISP(1:LBISP)/DBLE(NFRAMES)
      PRINT *,'Finished calculating the bispectrum.'
      RETURN
      END SUBROUTINE BISPECTRUM
C ******************************************************************************
      SUBROUTINE GENERATEBETA(M,N,MW,DPHI,DBETA)
C  Generate principal values of the arguments of the bispectrum of given phase
C  (beta) for testing the recursion algorithm.
C
      INTEGER, INTENT(IN) :: M,N,MW
      INTEGER :: I1,I2,J1,J2,K
      DOUBLE PRECISION, INTENT(IN) :: DPHI(0:M-1,0:N-1)
      DOUBLE PRECISION, INTENT(OUT) :: DBETA((MW+1)*(2*N-MW)*M*(M+2)/8)
C  Interface:
C  ==========
      INTERFACE
      FUNCTION BISPOS(I1,J1,I2,J2,M,N)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1,J1,I2,J2,M,N
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      END INTERFACE
      PRINT *,'GENERATE BETA'
      K=1
      DO J2=0,MW
        DO I2=0,M-1
          DO J1=0,N-1-J2
            DO I1=0,MIN(I2,M-1-I2)
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
      SUBROUTINE PHASERECURSION(M,N,MW,DBETA,DPHI)
C  Recursive algorithm to solve the phase equations.
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N,MW
      INTEGER :: I,J,K,I1,I2,J1,J2
      DOUBLE PRECISION :: R
      DOUBLE PRECISION, INTENT(OUT) :: DPHI(0:M-1,0:N-1)
      DOUBLE PRECISION, INTENT(IN) :: DBETA((MW+1)*(2*N-MW)*M*(M+2)/8)
C  Interface:
C  ==========
      INTERFACE
      FUNCTION BISPOS(I1,J1,I2,J2,M,N)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1,J1,I2,J2,M,N
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      END INTERFACE
C  Statements:
C  ===========
      DPHI=DBLE(0.0)
      DO K=2,M+N-2
        DO I=MAX(0,K+1-M),MIN(K,M-1)
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
     &            DBETA(BISPOS(I1,J1,I2,J2,M,N)))/R
              END IF
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE PHASERECURSION
C ******************************************************************************
      INTEGER FUNCTION BISPOS(I1,J1,I2,J2,M,N)
C  Calculate position in bispectrum array according to positions of the
C  frequency components in the spectrum matrices.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1,J1,I2,J2,M,N
      INTEGER :: K,L
      K=1+M*(M+2)*J2*(2*N-J2+1)/8
      IF (I2 .LT. M/2) THEN
        BISPOS=K+I2*(I2+1)*(N-J2)/2+J1*(I2+1)+I1
        RETURN
      END IF
      IF (I2 .EQ. M/2) THEN
        BISPOS=K+M*(M+2)*(N-J2)/8+J1*M/2+I1
        RETURN
      END IF
      IF (I2 .GT. M/2) THEN
        BISPOS=K+M*(M+2)*(N-J2)/8+
     &    (3*M-2*I2+2)*(2*I2-M)*(N-J2)/8+J1*(M-I2)+I1
        RETURN
      END IF
      END FUNCTION BISPOS
C *****************************************************************************
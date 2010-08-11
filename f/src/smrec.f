      SUBROUTINE BISPECTRUM(FILENAME,FPIXELS,LPIXELS,DAVG,DR,FTMETHOD,
     &  P,MW,BUFFERSIZE,ZBISP)
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
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3),P,MW,BUFFERSIZE
      INTEGER :: I1,I2,I3,J1,J2,K,LF,LT,LR,L,M,N,NPIXELS,NFRAMES,
     &  LBUFFER,INFO,LBISP
      INTEGER*8 :: PLAN
      DOUBLE PRECISION, INTENT(IN) :: DR,
     &  DAVG(LPIXELS(1)-FPIXELS(1)+1,LPIXELS(2)-FPIXELS(2)+1)
      DOUBLE PRECISION, INTENT(OUT) :: ZBISP(*)
      DOUBLE PRECISION, DIMENSION(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(2)-FPIXELS(2)+1) :: DBG,DIMG,DB
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:,:,:)
      DOUBLE COMPLEX, DIMENSION(P,P) :: ZIN,ZOUT,ZSP
      CHARACTER, INTENT(IN) :: FILENAME*256,FTMETHOD*10
C
C  Statements:
C  ===========
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=M*N
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
C  The padding size must be even:
      P=INT(CEILING(DBLE(P)/DBLE(2))*DBLE(2))
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
      NBISP=MW*(2*P-MW-1)*P*P/8
      ZBISP(1:NBISP)=CMPLX(0)
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
          ZIN(1:M,1:N)=CMPLX(DIMG)
          CALL ZIFFTSHIFT(P,P,ZIN)
          CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
          ZSP=ZOUT
          I3=1
          DO J2=1,MW
            DO I2=1,P-1
              DO J1=1,P-J2
                DO I1=1,MIN(I2,P-I2)
                  ZBISP(I3)=ZBISP(I3)+ZSP(I1,J1)*ZSP(I2,J2)*
     &              CONJG(ZSP(I1+I2,J1+J2))
                  PRINT *,I3
                  PRINT *,BISPOS(I1,J1,I2,J2,P)
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
C *****************************************************************************
      SUBROUTINE PHASERECURSION(P,MW,DBETA,DPHI)
C  Recursive algorithm to solve the phase equations.
C
C  Declarations:
C  =============
      INTEGER, INTENT(IN) :: P,MW
      DOUBLE PRECISION, INTENT(OUT) :: DPHI(P,P)
      DOUBLE PRECISION, INTENT(IN) :: DBETA(MW*(2*P-MW-1)*P*P/8)
      DPHI(1,1)=0
      DPHI(1,2)=0
      DPHI(2,1)=0
      
      RETURN
      END SUBROUTINE PHASERECURSION
C ******************************************************************************
      INTEGER FUNCTION BISPOS(I1,J1,I2,J2,P)
      INTEGER :: K,L
      K=0
      DO L=1,J2-1
        K=K+L*L
      END DO
      K=((2*P-1)*J2*(J2-1)-2*K)*P*P/16
      IF (I2 .LT. P/2+1) THEN
        K=K+I2*(I2-1)*(P-J2)/2+I2*(J1-1)+I1
      END IF
      IF (I2 .EQ. P/2+1) THEN
        K=K+P*(P+2)*(P-J2)/8+I2*(P-4)/2+I1
      END IF
      IF (I2 .GT. P/2+1) THEN
        K=K+(3*P-2*I2-2)*(2*I2-P-4)*(P-J2)/8+I2*(P-I2-1)+I1
      END IF
      BISPOS=K
      RETURN
      END FUNCTION BISPOS
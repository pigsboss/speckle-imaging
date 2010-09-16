      PROGRAM ISA
C  Iterative shift-and-add routine.
C  ================================
C
C  Usage:
C  ======
C  isa filename [-range=m_1,n_1] [-range=m_2,n_2] [...] -ref=filename
C    -init=filename -snr=SNR [-prefix=output]
C
C  Argument:
C  =========
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INTEGER, PARAMETER :: MAXNRNG=100,DEFAULTNUMIT=10
      INTEGER :: STATUS,UNIT,RNG(2,MAXNRNG),NAXES(3),NARGS,K,NPIXELS,
     &  L,NRNG,NUMIT
      DOUBLE PRECISION, PARAMETER :: DEFAULTSNR=1.0D6
      DOUBLE PRECISION :: DSNR
      CHARACTER(LEN=256) :: INFILE,PREFIX,ARG,BASENAME,EXTNAME,IFILE,
     &  RFILE
      INTERFACE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE ITERATIVESHIFTADD(INFILE,NRNG,RNG,IFILE,RFILE,DSNR,
     &  NUMIT,PREFIX)
      INTEGER, INTENT(IN) :: NRNG,RNG(2,NRNG),NUMIT
      DOUBLE PRECISION, INTENT(IN) :: DSNR
      CHARACTER(LEN=*), INTENT(IN) :: INFILE,IFILE,RFILE,PREFIX
      END SUBROUTINE ITERATIVESHIFTADD
      END INTERFACE
C  Statements:
C  ===========
      STATUS=0
C    Resolve the command line options:
C    =================================
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'isa filename [-range=m_1,n_1] [-range=m_2,n_2]'//
     &    ' [...] -ref=filename -init=filename [-snr=SNR]'//
     &    ' [-n=numit] [-prefix=output]'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,INFILE)
      CALL IMAGESIZE(INFILE,NAXES)
      RNG(:,1)=(/1,NAXES(3)/)
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)
      NRNG=0
      RFILE=''
      IFILE=''
      NUMIT=DEFAULTNUMIT
      DSNR=DEFAULTSNR
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF((INDEX(ARG,'-range=') .GT. 0).AND.(NRNG.LT.MAXNRNG))THEN
          NRNG=NRNG+1
          READ(ARG(INDEX(ARG,'-range=')+7:),*)
     &      RNG(1,NRNG),RNG(2,NRNG)
        ELSE IF((INDEX(ARG,'-range=') .GT. 0).AND.(NRNG.GE.MAXNRNG))THEN
          PRINT *,'Too many range definitions.'
          STOP
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=ARG(INDEX(ARG,'-prefix=')+8:)
        ELSE IF(INDEX(ARG,'-ref=').GT.0)THEN
          RFILE=ARG(INDEX(ARG,'-ref=')+5:)
        ELSE IF(INDEX(ARG,'-init=').GT.0)THEN
          IFILE=ARG(INDEX(ARG,'-init=')+6:)
        ELSE IF(INDEX(ARG,'-snr=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-snr=')+5:),*)DSNR
        ELSE IF(INDEX(ARG,'-n=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-n=')+3:),*)NUMIT
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          STOP
        END IF
      END DO
      IF(LEN_TRIM(IFILE) .LT. 1)THEN
        PRINT *,'the initial estimate must be specified.'
        STOP
      END IF
      IF(LEN_TRIM(RFILE) .LT. 1)THEN
        PRINT *,'the reference must be specified.'
        STOP
      END IF
      IF(NRNG .EQ. 0)THEN
        NRNG=1
      END IF
      WRITE(*,'(A,ES7.1)')' SNR=',DSNR
      WRITE(*,'(A,I3)')' number of iterations: ',NUMIT
      PRINT *,'input: '//TRIM(INFILE)
      PRINT *,'initial estimate: '//TRIM(IFILE)
      PRINT *,'reference: '//TRIM(RFILE)
      DO L=1,NRNG
        WRITE(*,'(A,I3,A,I5,A,I5)')' range ',L,': from ',RNG(1,L),
     &    ' to ',RNG(2,L)
      END DO
      PRINT *,'prefix of output: '//TRIM(PREFIX)
      CALL ITERATIVESHIFTADD(INFILE,NRNG,RNG,IFILE,RFILE,DSNR,NUMIT,
     &  PREFIX)
      STOP
      END PROGRAM ISA

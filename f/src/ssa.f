      PROGRAM SSA
C  Simple shift-and-add routine.
C  =============================
C
C  Usage:
C  ======
C  ssa filename [-range=m_1,n_1] [-range=m_2,n_2] [...] [-prefix=output]
C
C  Argument:
C  =========
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INTEGER, PARAMETER :: MAXNRNG=100
      INTEGER :: STATUS,UNIT,RNG(2,MAXNRNG),NAXES(3),NARGS,K,NPIXELS,
     &  L,NRNG
      CHARACTER(LEN=256) :: INFILE,PREFIX,ARG,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE SHIFTADD(INFILE,NRNG,RNG,PREFIX)
      INTEGER, INTENT(IN) :: NRNG,RNG(2,NRNG)
      CHARACTER(LEN=*), INTENT(IN) :: INFILE,PREFIX
      END SUBROUTINE SHIFTADD
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
        PRINT *,'ssa filename [-range=m_1,n_1] [-range=m_2,n_2]'//
     &      ' [...] [-prefix=output]'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,INFILE)
      CALL IMAGESIZE(INFILE,NAXES)
      RNG(:,1)=(/1,NAXES(3)/)
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)
      NRNG=0
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
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          STOP
        END IF
      END DO
      IF(NRNG .EQ. 0)THEN
        NRNG=1
      END IF
      PRINT *,'input: '//TRIM(INFILE)
      DO L=1,NRNG
        WRITE(*,'(A,I3,A,I5,A,I5)')' range ',L,': from ',RNG(1,L),
     &    ' to ',RNG(2,L)
      END DO
      PRINT *,'prefix of output: '//TRIM(PREFIX)
      CALL SHIFTADD(INFILE,NRNG,RNG,PREFIX)
      STOP
      END PROGRAM SSA

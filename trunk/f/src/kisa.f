      PROGRAM KISA
C  Kalman filter iterative shift-and-add routine.
C  ==============================================
C
C  Usage:
C  ======
C  kisa filename [-len-branch=n] [-num-branch=n] -psf=filename
C    [-init=filename] [-prefix=output] [-n=numit]
C
C  Argument:
C  =========
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INTEGER, PARAMETER :: DEFAULTNUMIT=10,DEFAULTLBR=1000,
     &  DEFAULTNBR=10
      INTEGER :: STATUS,UNIT,LBR,NBR,NAXES(3),NARGS,K,NPIXELS,
     &  L,NRNG,NUMIT
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
      SUBROUTINE KALMANITERATIVESHIFTADD(TARFILE,LBR,NBR,
     &  INIFILE,PSFFILE,NUMIT,PREFIX)
      INTEGER, INTENT(IN) :: LBR,NBR,NUMIT
      CHARACTER(LEN=*), INTENT(IN) :: TARFILE,INIFILE,PSFFILE,PREFIX
      END SUBROUTINE KALMANITERATIVESHIFTADD
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
        PRINT *,'kisa filename [-len-branch=n] [-num-branch=n]'//
     &    ' -psf=filename [-init=filename] [-prefix=output] [-n=numit]'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,INFILE)
      CALL IMAGESIZE(INFILE,NAXES)
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)//'_kisa'
      LBR=DEFAULTLBR
      NBR=DEFAULTNBR
      RFILE=''
      IFILE=''
      NUMIT=DEFAULTNUMIT
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-len-branch=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-len-branch=')+12:),*)LBR
        ELSE IF(INDEX(ARG,'-num-branch=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-num-branch=')+12:),*)NBR
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=ARG(INDEX(ARG,'-prefix=')+8:)
        ELSE IF(INDEX(ARG,'-psf=').GT.0)THEN
          RFILE=ARG(INDEX(ARG,'-psf=')+5:)
        ELSE IF(INDEX(ARG,'-init=').GT.0)THEN
          IFILE=ARG(INDEX(ARG,'-init=')+6:)
        ELSE IF(INDEX(ARG,'-n=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-n=')+3:),*)NUMIT
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          STOP
        END IF
      END DO
      IF(LEN_TRIM(IFILE) .LT. 1)THEN
        PRINT *,'the reference has not been specified.'
      END IF
      IF(LEN_TRIM(RFILE) .LT. 1)THEN
        PRINT *,'the reference must be specified.'
        STOP
      END IF
      WRITE(*,'(A,I2)')' number of iterations: ',NUMIT
      PRINT *,'input: '//TRIM(INFILE)
      PRINT *,'reference: '//TRIM(RFILE)
      PRINT *,'prefix of output: '//TRIM(PREFIX)
      CALL KALMANITERATIVESHIFTADD(INFILE,LBR,NBR,
     &  IFILE,RFILE,NUMIT,PREFIX)
      STOP
      END PROGRAM KISA
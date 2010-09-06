      PROGRAM SM
C  sm filename_target filename_ref -level=y_2_max
C    [-target-range=m_1,n_1] [-target-range=m_2,n_2] [...]
C    [-ref-range=m_1,n_1] [-ref-range=m_1,n_1] [...]
C    [-psd=filename] [-prefix=output]
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: MAXNRNG=100
      INTEGER :: Y2MAX,NAXES(3),NARGS,NTRNG,NRRNG,TRNG(2,MAXNRNG),
     &  RRNG(2,MAXNRNG),K
      CHARACTER(LEN=256) :: TFILE,RFILE,PREFIX,PSDFILE,
     &  ARG,BASENAME,EXTNAME
C
      INTERFACE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE SPECKLEMASKING(TFILE,NTRNG,TRNG,RFILE,NRRNG,RRNG,
     &  Y2MAX,PSDFILE,PREFIX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NTRNG,TRNG(2,NTRNG),NRRNG,RRNG(2,NRRNG),
     &  Y2MAX
      CHARACTER(LEN=*), INTENT(IN) :: TFILE,RFILE,PSDFILE,PREFIX
      END SUBROUTINE SPECKLEMASKING
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      END INTERFACE
C
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'sm filename_target filename_ref -level=y_2_max'//
     &    ' [-target-range=m_1,n_1] [-target-range=m_2,n_2] [...]'//
     &    ' [-ref-range=m_1,n_1] [-ref-range=m_1,n_1] [...]'//
     &    ' [-prefix=output]'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,TFILE)
      CALL RESOLVEPATH(TFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)//'_sm'
      CALL IMAGESIZE(TFILE,NAXES)
      TRNG(:,1)=(/1,NAXES(3)/)
      CALL GET_COMMAND_ARGUMENT(2,RFILE)
      CALL IMAGESIZE(RFILE,NAXES)
      RRNG(:,1)=(/1,NAXES(3)/)
      NTRNG=0
      NRRNG=0
      Y2MAX=0
      PSDFILE=''
      DO K=3,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF((INDEX(ARG,'-target-range=') .GT. 0).AND.
     &    (NTRNG.LT.MAXNRNG))THEN
          NTRNG=NTRNG+1
          READ(ARG(INDEX(ARG,'-target-range=')+14:),*)
     &      TRNG(1,NTRNG),TRNG(2,NTRNG)
        ELSE IF((INDEX(ARG,'-target-range=') .GT. 0).AND.
     &    (NTRNG.GE.MAXNRNG))THEN
          PRINT *,'error: too many target range definitions.'
          STOP
        ELSE IF((INDEX(ARG,'-ref-range=') .GT. 0).AND.
     &    (NRRNG.LT.MAXNRNG))THEN
          NRRNG=NRRNG+1
          READ(ARG(INDEX(ARG,'-ref-range=')+11:),*)
     &      RRNG(1,NRRNG),RRNG(2,NRRNG)
        ELSE IF((INDEX(ARG,'-ref-range=') .GT. 0).AND.
     &    (NRRNG.GE.MAXNRNG))THEN
          PRINT *,'error: too many reference range definitions.'
          STOP
        ELSE IF(INDEX(ARG,'-level=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-level=')+7:),*)Y2MAX
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=ARG(INDEX(ARG,'-prefix=')+8:)
        ELSE IF(INDEX(ARG,'-psd=').GT.0)THEN
          PSDFILE=ARG(INDEX(ARG,'-psd=')+5:)
        ELSE
          PRINT *,'error: unknown argument '//TRIM(ARG)
          STOP
        END IF
      END DO
      IF(Y2MAX .EQ. 0)THEN
        PRINT *,'error: bispectral level invalid.'
        STOP
      END IF
      IF(NTRNG .EQ. 0)THEN
        NTRNG=1
      END IF
      IF(NRRNG .EQ. 0)THEN
        NRRNG=1
      END IF
      CALL SPECKLEMASKING(TFILE,NTRNG,TRNG,RFILE,NRRNG,RRNG,
     &  Y2MAX,PSDFILE,PREFIX)
      STOP
      END PROGRAM SM

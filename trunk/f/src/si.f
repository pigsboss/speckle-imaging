      PROGRAM SI
C  Speckle interferometry.
C  =======================
C
C  Usage:
C  ======
C  si filename_target filename_reference [-target-range=m,n] [-ref-range=m,n]
C    [-prefix=...]
C
C  Arguments:
C  ==========
C
      IMPLICIT NONE
      INTEGER :: NAXES(3),TRNG(2),RRNG(2),NARGS,K,PX,PY
      CHARACTER(LEN=256) :: TFILE,RFILE,PREFIX,ARG,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      INTEGER, INTENT(OUT) :: NAXES(3)
      END SUBROUTINE IMAGESIZE
C
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER*(*), INTENT(IN) :: PATH
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
C
      SUBROUTINE SPECKLEINTERFEROMETRY(TFILE,TRNG,RFILE,RRNG,PX,PY,
     &  PREFIX,BUFFERSIZE)
      INTEGER, INTENT(IN) :: TRNG(2),RRNG(2),PX,PY
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      CHARACTER*(*), INTENT(IN) :: TFILE,RFILE,PREFIX
      END SUBROUTINE SPECKLEINTERFEROMETRY
      END INTERFACE
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'si filename_target filename_reference '//
     &    '[-target-range=m,n] [-ref-range=m,n] [-prefix=...]'
        STOP
      END IF
      PX=0
      PY=0
      CALL GET_COMMAND_ARGUMENT(1,TFILE)
      CALL GET_COMMAND_ARGUMENT(2,RFILE)
      CALL IMAGESIZE(TFILE,NAXES)
      PX=MAX(PX,NAXES(1))
      PY=MAX(PY,NAXES(2))
      TRNG=(/1,NAXES(3)/)
      CALL IMAGESIZE(RFILE,NAXES)
      PX=MAX(PX,NAXES(1))
      PY=MAX(PY,NAXES(2))
      RRNG=(/1,NAXES(3)/)
      PX=INT(IDNINT(
     &  DEXP(DLOG(2.0D0)*CEILING(DLOG(DBLE(PX))/DLOG(2.0D0)))))
      PY=INT(IDNINT(
     &  DEXP(DLOG(2.0D0)*CEILING(DLOG(DBLE(PY))/DLOG(2.0D0)))))
      CALL RESOLVEPATH(TFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)
      CALL RESOLVEPATH(RFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(PREFIX)//'_'//TRIM(BASENAME)//'_si'
      DO K=3,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-target-range=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-target-range=')+14:),*) TRNG(1),TRNG(2)
        ELSE IF(INDEX(ARG,'-ref-range=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-ref-range=')+11:),*) RRNG(1),RRNG(2)
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=ARG(INDEX(ARG,'-prefix=')+8:)
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          RETURN
        END IF
      END DO
      CALL SPECKLEINTERFEROMETRY(TFILE,TRNG,RFILE,RRNG,PX,PY,PREFIX)
      STOP
      END PROGRAM SI

      INCLUDE './lib/math.f'
      INCLUDE './lib/reconstruct.f'
      INCLUDE './lib/image.f'
      INCLUDE './lib/signal.f'
C ******************************************************************************
      SUBROUTINE PRINTERROR(STATUS)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: STATUS
      CHARACTER :: ERRTEXT*30,ERRMESSAGE*80
      IF (STATUS .LE. 0)RETURN
      CALL FTGERR(STATUS,ERRTEXT)
      PRINT *,'FITSIO ERROR STATUS =',STATUS,': ',ERRTEXT
      CALL FTGMSG(ERRMESSAGE)
      DO WHILE (ERRMESSAGE .NE. ' ')
        PRINT *,ERRMESSAGE
        CALL FTGMSG(ERRMESSAGE)
      ENDDO
      RETURN
      END SUBROUTINE PRINTERROR
C ******************************************************************************
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: STATUS
      INTEGER :: UNIT,BLOCKSIZE
      CHARACTER*(*), INTENT(IN) :: FILENAME
      INTERFACE
      SUBROUTINE PRINTERROR(STATUS)
      INTEGER, INTENT(IN) :: STATUS
      END SUBROUTINE PRINTERROR
      END INTERFACE
      IF (STATUS .GT. 0)RETURN
      CALL FTGIOU(UNIT,STATUS)
      CALL FTOPEN(UNIT,FILENAME,1,BLOCKSIZE,STATUS)
      IF (STATUS .EQ. 0)THEN
        CALL FTDELT(UNIT,STATUS)
      ELSE IF (STATUS .EQ. 103)THEN
        STATUS=0
        CALL FTCMSG
      ELSE
        STATUS=0
        CALL FTCMSG
        CALL FTDELT(UNIT,STATUS)
      END IF
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE DELETEFILE
C ******************************************************************************
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(IN) :: PATH
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      CHARACTER :: FILENAME*80
      INTEGER :: K,L
      K=SCAN(PATH,'/',.TRUE.)
      L=LEN_TRIM(PATH)
      READ(PATH(K+1:L),*) FILENAME
      K=SCAN(FILENAME,'.')
      L=LEN_TRIM(FILENAME)
      READ(FILENAME(1:K-1),*) BASENAME
      K=SCAN(FILENAME,'.',.TRUE.)
      READ(FILENAME(K+1:L),*) EXTNAME
      RETURN
      END SUBROUTINE RESOLVEPATH
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************

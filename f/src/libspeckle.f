      SUBROUTINE PRINTERROR(STATUS)
      INTEGER :: STATUS
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
      ENDSUBROUTINE PRINTERROR
C ******************************************************************************
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      INTEGER :: STATUS,UNIT,BLOCKSIZE
      CHARACTER*(*) :: FILENAME
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
      RETURN
      ENDSUBROUTINE DELETEFILE
C ******************************************************************************
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER :: STATUS,UNIT,BLOCKSIZE,BITPIX,NAXIS,GROUP
      INTEGER :: FPIXELS(*),LPIXELS(*),NAXES(3)
      DOUBLE PRECISION :: ARRAY(*)
      CHARACTER*(*) :: FILENAME
      STATUS=0
      CALL DELETEFILE(FILENAME,STATUS)
      CALL FTGIOU(UNIT,STATUS)
      BLOCKSIZE=1
      BITPIX=-64
      NAXIS=3
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      CALL FTINIT(UNIT,FILENAME,BLOCKSIZE,STATUS)
      CALL FTPHPS(UNIT,BITPIX,NAXIS,NAXES,STATUS)
      GROUP=1
      CALL FTPSSD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,ARRAY,STATUS)
      CALL FTCLOS(UNIT,STATUS)
      CALL FTFIOU(UNIT,STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      ENDSUBROUTINE WRITEIMAGE
C ******************************************************************************
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER :: STATUS,UNIT,READWRITE,BLOCKSIZE,NAXIS,NFOUND
      INTEGER :: GROUP,FPIXELS(*),LPIXELS(*),INCS(3),NAXES(3)
      DOUBLE PRECISION :: ARRAY(*),NULLVAL
      LOGICAL :: ANYF
      CHARACTER*(*) :: FILENAME
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      READWRITE=0
      CALL FTOPEN(UNIT,FILENAME,READWRITE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'READIMAGE FAILED TO READ THE NAXISN KEYWORDS.'
        RETURN
      ENDIF
      GROUP=1
      NAXIS=3
      INCS(1)=1
      INCS(2)=1
      INCS(3)=1
      NULLVAL=-999
      CALL FTGSVD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,INCS,
     &  NULLVAL,ARRAY,ANYF,STATUS)
      IF (ANYF)THEN
        PRINT *,'One or more pixels are undefined in the FITS image.'
        RETURN
      ENDIF
      CALL FTCLOS(UNIT, STATUS)
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      ENDSUBROUTINE READIMAGE
C ******************************************************************************
      SUBROUTINE AVERAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER :: FPIXELS(*),LPIXELS(*),K,NPIXELS,F(3),L(3)
      CHARACTER*(*) :: FILENAME
      DOUBLE PRECISION :: ARRAY(*),DA
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:)
      INTEGER :: NAXES(3)
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      F(1)=FPIXELS(1)
      F(2)=FPIXELS(2)
      L(1)=LPIXELS(1)
      L(2)=LPIXELS(2)
      ALLOCATE(BUFFER(NPIXELS))
      DO K=1,NPIXELS
        ARRAY(K)=0
      ENDDO
      DA=1.0
      DO K=FPIXELS(3),LPIXELS(3)
        F(3)=K
        L(3)=K
        CALL READIMAGE(FILENAME,F,L,BUFFER)
        CALL DAXPY(NPIXELS,DA,BUFFER,1,ARRAY,1)
      ENDDO
      DA=1.0/DBLE(NAXES(3))
      CALL DSCAL(NPIXELS,DA,ARRAY,1)
      DEALLOCATE(BUFFER)
      RETURN
      ENDSUBROUTINE AVERAGE
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
C ******************************************************************************
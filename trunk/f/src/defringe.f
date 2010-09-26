      PROGRAM DEFRINGE
C  Usage:
C  ======
C  defringe filename_img [-from=n] [-to=n] [-width=w]
C    [-o=output] [-buffer=n]
C
C  Purpose:
C  ========
C  cut off the fringe level of each frame.
C
      IMPLICIT NONE
      INTEGER,PARAMETER :: DEFAULTWIDTH=2,DEFAULTBUFSIZ=128
      INTEGER :: STATUS,NARGS,K,X,Y,NAXES(3),BUFSIZ,LBUF,NBUF,WIDTH,
     &  NFRMS,L,L1,L2,FROM,TO
      DOUBLE PRECISION :: DLEV
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:),DFLAG(:,:)
      CHARACTER(LEN=256) :: INFILE,OUTFILE,ARG,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER*(*), INTENT(IN) :: PATH
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      INTEGER, INTENT(INOUT) :: STATUS
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE DELETEFILE
      SUBROUTINE READIMAGE(IMGFILE,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      END SUBROUTINE READIMAGE
      SUBROUTINE APPENDIMAGE(IMGFILE,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      END SUBROUTINE APPENDIMAGE
      SUBROUTINE WRITEIMAGE(IMGFILE,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      END SUBROUTINE WRITEIMAGE
      SUBROUTINE PRINTERROR(STATUS)
      INTEGER, INTENT(IN) :: STATUS
      END SUBROUTINE PRINTERROR
      END INTERFACE
      STATUS=0
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'defringe filename [-o=filename] [-width=w] [-from=n]'
        PRINT *,'  [-to=n] [-buffer=buffer_size]'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,INFILE)
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
      OUTFILE=TRIM(BASENAME)//'_def.fits'
      CALL IMAGESIZE(INFILE,NAXES)
      BUFSIZ=DEFAULTBUFSIZ
      WIDTH=DEFAULTWIDTH
      FROM=1
      TO=NAXES(3)
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-o=').GT.0)THEN
          OUTFILE=ARG(INDEX(ARG,'-o=')+3:)
        ELSE IF(INDEX(ARG,'-from=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-from=')+6:),*) FROM
        ELSE IF(INDEX(ARG,'-to=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-to=')+4:),*) TO
        ELSE IF(INDEX(ARG,'-buffer=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-buffer=')+8:),*) BUFSIZ
        ELSE IF(INDEX(ARG,'-width=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-width=')+7:),*) WIDTH
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          RETURN
        END IF
      END DO
      NFRMS=TO-FROM+1
      LBUF=INT(FLOOR(DBLE(BUFSIZ*1024*1024)/DBLE(8*NAXES(1)*NAXES(2))))
      NBUF=INT(CEILING(DBLE(NFRMS)/DBLE(LBUF)))
      WRITE(*,*) 'input: '//TRIM(INFILE)
      WRITE(*,*) 'output: '//TRIM(OUTFILE)
      WRITE(*,'(A,I5,A,I5)') ' from frame ',FROM,' to frame ',TO
      WRITE(*,'(A,I3,A)') ' buffer size: ',BUFSIZ,'MB'
      WRITE(*,'(A,I3,A)') ' fringe width: ',WIDTH,' pixels'
C
      CALL DELETEFILE(OUTFILE,STATUS)
      ALLOCATE(DBUF(NAXES(1),NAXES(2),LBUF),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        WRITE(*,*)'error: out of memory.'
        STOP
      END IF
      ALLOCATE(DFLAG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        STOP
      END IF
C
      DO X=1,NAXES(1)
        DO Y=1,NAXES(2)
          IF((X .LE. WIDTH) .OR. (X .GE. NAXES(1)+1-WIDTH) .OR.
     &      (Y .LE. WIDTH) .OR. (Y .GE. NAXES(2)+1-WIDTH))THEN
            DFLAG(X,Y)=1.0D0
          ELSE
            DFLAG(X,Y)=0.0D0
          END IF
        END DO
      END DO
      DO K=1,NBUF
        L1=(K-1)*LBUF+FROM
        L2=MIN(K*LBUF,TO)
        CALL READIMAGE(INFILE,(/1,1,L1/),(/NAXES(1),NAXES(2),L2/),DBUF)
        DO L=1,L2-L1+1
          DLEV=MAXVAL(DFLAG*DBUF(:,:,L))
          WRITE(*,'(A,I5,A,ES7.1)')' frame ',L1+L-1,
     &      ', fringe level: ',DLEV
          DBUF(:,:,L)=DBUF(:,:,L)-DLEV
          DO X=1,NAXES(1)
            DO Y=1,NAXES(2)
              IF(DBUF(X,Y,L) .LT. 0.0D0)THEN
                DBUF(X,Y,L)=0.0D0
              END IF
            END DO
          END DO
        END DO
        CALL APPENDIMAGE(OUTFILE,(/1,1,L1/),
     &    (/NAXES(1),NAXES(2),L2/),DBUF)
      END DO
      DEALLOCATE(DBUF)
      DEALLOCATE(DFLAG)
      STOP
      END PROGRAM DEFRINGE

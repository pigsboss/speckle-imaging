      PROGRAM BGCUT
C  Usage:
C  ======
C  bgcut filename_img filename_bg [-from=] [-to=] [-r=radius] [-prefix=]
C    [-buffer=]
C
C  Purpose:
C  ========
C  cut off the previously calculated background from the images.
C
C  Arguments:
C  ==========
C  filename_img - multi-frame image filename (flat field corrected).
C  filename_bg  - background filename (previously calculated).
C  radius       - radius of the signal region. 
c                   [default=0.5*MIN(NAXES(1),NAXES(2))]
C  from         - indicates the first frame. [default=1]
C  to           - indicates the last frame. [default=NAXES(3)]
C  prefix       - prefix of output filenames in this program.
C                   [default: basename of filename_img]
C  buffer       - buffer size in MB. [default=32]
C
      IMPLICIT NONE
      INTEGER :: NAXES(3),K,X,Y,NARGS,FROM,TO,NFRAMES,BUFFERSIZE,
     &  LBUFFER,NBUFFER,L1,L2,L,STATUS
      DOUBLE PRECISION :: DR
      DOUBLE PRECISION :: XC,YC,DX,DY
      DOUBLE PRECISION, ALLOCATABLE :: DBG(:,:),DSTD(:,:),
     &  DFLAG(:,:),DBUF(:,:,:)
      CHARACTER(LEN=256) :: IMGFILE,BGFILE,ARG,PREFIX,BASENAME,EXTNAME
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
      CALL GET_COMMAND_ARGUMENT(1,IMGFILE)
      CALL GET_COMMAND_ARGUMENT(2,BGFILE)
      CALL RESOLVEPATH(IMGFILE,BASENAME,EXTNAME)
      PREFIX=BASENAME
      CALL IMAGESIZE(IMGFILE,NAXES)
      BUFFERSIZE=32
      FROM=1
      TO=NAXES(3)
      DR=0.5*DBLE(MIN(NAXES(1),NAXES(2)))
      DO K=3,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-r=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-r=')+3:256),*) DR
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-prefix=')+8:256),*) PREFIX
        ELSE IF(INDEX(ARG,'-from=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-from=')+6:256),*) FROM
        ELSE IF(INDEX(ARG,'-to=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-to=')+4:256),*) TO
        ELSE IF(INDEX(ARG,'-buffer=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-buffer=')+8:256),*) BUFFERSIZE
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          RETURN
        END IF
      END DO
      NFRAMES=TO-FROM+1
      LBUFFER=INT(FLOOR(DBLE(BUFFERSIZE*1024*1024)/
     &  DBLE(8*NAXES(1)*NAXES(2))))
      NBUFFER=INT(CEILING(DBLE(NFRAMES)/DBLE(LBUFFER)))
      WRITE(*,'(A,F6.1)') ' radius: ',DR
      WRITE(ARG,*) FROM
      WRITE(*,'(A,A)') ' from: ',TRIM(ADJUSTL(ARG))//' frame'
      WRITE(ARG,*) TO
      WRITE(*,'(A,A)') ' to: ',TRIM(ADJUSTL(ARG))//' frame'
      WRITE(ARG,*) BUFFERSIZE
      WRITE(*,'(A,A)') ' buffer size: ',TRIM(ADJUSTL(ARG))//'MB'
      WRITE(*,'(A,A)') ' output filename prefix: ',TRIM(PREFIX)
      CALL DELETEFILE(TRIM(PREFIX)//'_sig.fits',STATUS)
      ALLOCATE(DBG(NAXES(2),NAXES(1)))
      ALLOCATE(DSTD(NAXES(2),NAXES(1)))
      ALLOCATE(DFLAG(NAXES(2),NAXES(1)))
      ALLOCATE(DBUF(NAXES(2),NAXES(1),LBUFFER))
      XC=0.5*DBLE(NAXES(1)+1)
      YC=0.5*DBLE(NAXES(2)+1)
      DO X=1,NAXES(1)
        DO Y=1,NAXES(2)
          DX=DBLE(X)-XC
          DY=DBLE(Y)-YC
          IF(DX*DX+DY*DY .GE. DR*DR)THEN
            DFLAG(Y,X)=DBLE(1)
          ELSE
            DFLAG(Y,X)=DBLE(0)
          END IF
        END DO
      END DO
      CALL READIMAGE(BGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
      DSTD=DBLE(0)
      DO K=1,NBUFFER
        L1=(K-1)*LBUFFER+FROM
        L2=MIN(K*LBUFFER,TO)
        CALL READIMAGE(IMGFILE,(/1,1,L1/),(/NAXES(1),NAXES(2),L2/),DBUF)
        DO L=1,L2-L1+1
          DBUF(1:NAXES(2),1:NAXES(1),L)=DBUF(1:NAXES(2),1:NAXES(1),L)-
     &      DBG
          DSTD=DSTD+DBUF(1:NAXES(2),1:NAXES(1),L)*
     &      DBUF(1:NAXES(2),1:NAXES(1),L)*DFLAG
          DO X=1,NAXES(1)
            DO Y=1,NAXES(2)
              IF(DBUF(Y,X,L).LT.DBLE(0))THEN
                DBUF(Y,X,L)=DBLE(0)
              END IF
            END DO
          END DO
        END DO
        CALL APPENDIMAGE(TRIM(PREFIX)//'_sig.fits',(/1,1,L1/),
     &    (/NAXES(1),NAXES(2),L2/),DBUF)
      END DO
      PRINT *,'signal file: ',TRIM(PREFIX)//'_sig.fits'
      DSTD=DSTD/DBLE(NFRAMES)
      WRITE(*,'(A,ES10.3)') ' square root of mean variance: ',
     &  DSQRT(SUM(DSTD)/SUM(DFLAG))
      DSTD=DSQRT(DSTD)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_dev.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DSTD)
      PRINT *,'deviation map: ',TRIM(PREFIX)//'_dev.fits'
      DEALLOCATE(DBUF)
      DEALLOCATE(DBG)
      DEALLOCATE(DSTD)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      STOP
      END PROGRAM BGCUT

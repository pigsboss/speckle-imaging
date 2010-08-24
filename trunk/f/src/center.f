      PROGRAM CENTER
C  Arguments:
C  ==========
C  NARGS  - Number of command arguments.
C  NAXES  - Lengths of axes of image.
C  STATUS - Status.
C  K      - Index of command arguments.
C  DIMG   - Input image to be centred and output.
C  INFILE - Input filename.
C  PREFIX - Prefix of output filename.
C  ARG    - Command argument.
      IMPLICIT NONE
      INTEGER :: NARGS,NAXES(3),STATUS,K
      DOUBLE PRECISION, ALLOCATABLE :: DIMG(:,:)
      CHARACTER(LEN=256) :: INFILE,PREFIX,ARG,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(IN) :: NAXES(3)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE CENTERIMAGE(NX,NY,DIMG)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(INOUT) :: DIMG(NX,NY)
      END SUBROUTINE CENTERIMAGE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
      END INTERFACE
      STATUS=0
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,INFILE)
      PRINT *,'input: ',TRIM(INFILE)
      CALL IMAGESIZE(INFILE,NAXES)
      WRITE(*,'(A,I4,A,I4)')' size (width x height):',NAXES(1),' x',
     &  NAXES(2)
      IF(NAXES(3).GT.1)THEN
        PRINT *,'input multi-frame FITS file.'
        PRINT *,'only process the first frame.'
      END IF
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=TRIM(ARG(INDEX(ARG,'-prefix=')+8:))
        ELSE
          PRINT *,'unknown argument '//ARG
          STOP
        END IF
      END DO
      ALLOCATE(DIMG(NAXES(2),NAXES(1)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(INFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      CALL CENTERIMAGE(NAXES(1),NAXES(2),DIMG)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_ctr.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DIMG)
      PRINT *,'output: '//TRIM(PREFIX)//'_ctr.fits'
      STOP
      END PROGRAM CENTER
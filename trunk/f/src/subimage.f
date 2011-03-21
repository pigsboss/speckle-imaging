      PROGRAM MAIN
C  Extract subimages.
C
C  Usage:
C  ======
C  subimage filename_img [-s] [-fpixels=n,n,n] [-lpixels=n,n,n] [-o=filename_sub]
C
C  Purpose:
C  ========
C  take filename_img as input. extract subimages defined by fpixels and lpixels
C  and write the output to file named filename_sub.
C
C  Arguments:
C  ==========
C  -s           - extract and split subimages.
C  filename_img - input
C  n,n,n        - integers
C  filename_sub - output
C
      IMPLICIT NONE
      INTEGER :: FPIXELS(3),LPIXELS(3),NAXES(3),K,NARGS
      CHARACTER(LEN=256) :: IMGFILE,SUBFILE,ARG,BASENAME,EXTNAME
      LOGICAL :: SPLIT
      INTERFACE
      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      END SUBROUTINE IMAGESIZE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER*(*), INTENT(IN) :: PATH
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE SUBIMAGE(IMGFILE,FPIXELS,LPIXELS,SUBFILE,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      CHARACTER*(*), INTENT(IN) :: IMGFILE,SUBFILE
      END SUBROUTINE SUBIMAGE
      SUBROUTINE SPLITIMAGE(IMGFILE,FPIXELS,LPIXELS,SUBFILE,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      CHARACTER*(*), INTENT(IN) :: IMGFILE,SUBFILE
      END SUBROUTINE SPLITIMAGE
      END INTERFACE
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,IMGFILE)
      CALL IMAGESIZE(IMGFILE,NAXES)
      FPIXELS=(/1,1,1/)
      LPIXELS=NAXES
      CALL RESOLVEPATH(IMGFILE,BASENAME,EXTNAME)
      SUBFILE=TRIM(BASENAME)//'_sub.fits'
      SPLIT=.FALSE.
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-fpixels=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-fpixels=')+9:),*)
     &      FPIXELS(1),FPIXELS(2),FPIXELS(3)
        ELSE IF(INDEX(ARG,'-lpixels=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-lpixels=')+9:),*)
     &      LPIXELS(1),LPIXELS(2),LPIXELS(3)
        ELSE IF(INDEX(ARG,'-o=').GT.0)THEN
          SUBFILE=ARG(INDEX(ARG,'-o=')+3:)
        ELSE IF(INDEX(ARG,'-s').GT.0)THEN
          SPLIT=.TRUE.
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          STOP
        END IF
      END DO
      WRITE(*,'(A)') ' input: '//TRIM(IMGFILE)
      WRITE(*,'(A)') ' output: '//TRIM(SUBFILE)
      WRITE(*,'(A,I3,A,I3)')
     &  ' subimage size (width x height): ',
     &  LPIXELS(1)-FPIXELS(1)+1,' x',LPIXELS(2)-FPIXELS(2)+1
      WRITE(*,'(A,I5,A,I5)') ' from ',FPIXELS(3),
     &  ' to ',LPIXELS(3)
      IF(SPLIT)THEN
        PRINT *,'Extract subimages into split files.'
        CALL SPLITIMAGE(IMGFILE,FPIXELS,LPIXELS,SUBFILE)
      ELSE
        PRINT *,'Extract subimages into a single file.'
        CALL SUBIMAGE(IMGFILE,FPIXELS,LPIXELS,SUBFILE)
      END IF
      STOP
      END PROGRAM MAIN

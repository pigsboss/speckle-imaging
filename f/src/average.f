      PROGRAM MAIN
C  Average images.
C
C  Usage:
C  ======
C  average filename_img [-fpixels=n,n,n] [-lpixels=n,n,n] [-o=filename_avg]
C
C  Purpose:
C  ========
C  take filename_img as input. average subimages defined by fpixels and lpixels
C  and write the output to file named filename_avg.
C
C  Arguments:
C  ==========
C  filename_img - input
C  n,n,n        - integers
C  filename_avg - output
C
      IMPLICIT NONE
      INTEGER :: FPIXELS(3),LPIXELS(3),NAXES(3),K,NARGS
      DOUBLE PRECISION, ALLOCATABLE :: DAVG(:,:)
      CHARACTER(LEN=256) :: FILENAME,ARG,AVGFILE,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE AVERAGE(FILENAME,FPIXELS,LPIXELS,ARRAY,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      DOUBLE PRECISION, INTENT(OUT) :: ARRAY(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(2)-FPIXELS(2)+1)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE AVERAGE
      SUBROUTINE WRITEIMAGE(IMGFILE,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: IMGFILE
      END SUBROUTINE WRITEIMAGE
      END INTERFACE
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,FILENAME)
      CALL IMAGESIZE(FILENAME,NAXES)
      FPIXELS=(/1,1,1/)
      LPIXELS=NAXES
      CALL RESOLVEPATH(FILENAME,BASENAME,EXTNAME)
      AVGFILE=TRIM(BASENAME)//'_avg.fits'
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-fpixels=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-fpixels=')+9:),*)
     &      FPIXELS(1),FPIXELS(2),FPIXELS(3)
        ELSE IF(INDEX(ARG,'-lpixels=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-lpixels=')+9:),*)
     &      LPIXELS(1),LPIXELS(2),LPIXELS(3)
        ELSE IF(INDEX(ARG,'-o=').GT.0)THEN
          AVGFILE=ARG(INDEX(ARG,'-o=')+3:)
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          RETURN
        END IF
      END DO
      WRITE(*,'(A)') ' input: '//TRIM(FILENAME)
      WRITE(*,'(A)') ' output: '//TRIM(AVGFILE)
      WRITE(*,'(A,I3,A,I3)')
     &  ' average subimage size (width x height): ',
     &  LPIXELS(1)-FPIXELS(1)+1,' x',LPIXELS(2)-FPIXELS(2)+1
      WRITE(*,'(A,I5,A,I5)') ' average from ',FPIXELS(3),
     &  ' to ',LPIXELS(3)
      ALLOCATE(DAVG(LPIXELS(1)+1-FPIXELS(1),LPIXELS(2)+1-FPIXELS(2)))
      CALL AVERAGE(TRIM(FILENAME),FPIXELS,LPIXELS,DAVG)
      CALL WRITEIMAGE(TRIM(AVGFILE),(/1,1,1/),
     &  (/LPIXELS(1)+1-FPIXELS(1),LPIXELS(2)+1-FPIXELS(2),1/),DAVG)
      DEALLOCATE(DAVG)
      STOP
      END PROGRAM

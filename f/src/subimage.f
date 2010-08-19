      PROGRAM MAIN
C  Extract subimages.
C
C  Usage:
C  ======
C  subimage filename_img [-fpixels=n,n,n] [-lpixels=n,n,n] [-o=filename_sub]
C
C  Purpose:
C  ========
C  take filename_img as input. extract subimages defined by fpixels and lpixels
C  and write the output to file named filename_sub.
C
C  Arguments:
C  ==========
C  filename_img - input
C  n,n,n        - integers
C  filename_sub - output
C
      IMPLICIT NONE
      INTEGER :: FPIXELS(3),LPIXELS(3),NAXES(3),K,NARGS
      CHARACTER(LEN=256) :: IMGFILE,SUBFILE,ARG,BASENAME,EXTNAME
      INTERFACE
      INCLUDE 'ispeckle.f'
      END INTERFACE
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,IMGFILE)
      CALL IMAGESIZE(IMGFILE,NAXES)
      FPIXELS=(/1,1,1/)
      LPIXELS=NAXES
      CALL RESOLVEPATH(IMGFILE,BASENAME,EXTNAME)
      SUBFILE=TRIM(BASENAME)//'_sub.fits'
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
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          RETURN
        END IF
      END DO
      CALL SUBIMAGE(IMGFILE,FPIXELS,LPIXELS,SUBFILE)
      STOP
      END PROGRAM MAIN

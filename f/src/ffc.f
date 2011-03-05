      PROGRAM FFC
C  Usage:
C  ======
C  ffc filename_target -flat=filename_flat -dark=filename_dark
C    [-fpixels=n,n,n] [-lpixels=n,n,n] [-o=filename_corrected]
C
C  Purpose:
C  ========
C  Flat field correction
C
      IMPLICIT NONE
      INTEGER :: FPIXELS(3),LPIXELS(3),NARGS,K,NAXES(3)
      CHARACTER(LEN=256) :: ARG,TARGFILE,FLATFILE,DARKFILE,CRCTFILE,
     &  BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE FLATFIELDCORRECT(TARGFILE,FLATFILE,DARKFILE,
     &  FPIXELS,LPIXELS,CRCTFILE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      CHARACTER*(*), INTENT(IN) :: TARGFILE,FLATFILE,DARKFILE,CRCTFILE
      END SUBROUTINE FLATFIELDCORRECT
      SUBROUTINE RESOLVEPATH(FILENAME,BASENAME,EXTNAME)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE IMAGESIZE
      END INTERFACE
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,TARGFILE)
      CALL RESOLVEPATH(TARGFILE,BASENAME,EXTNAME)
      CRCTFILE=TRIM(BASENAME)//'_crct.fits'
      CALL IMAGESIZE(TARGFILE,NAXES)
      FPIXELS=(/1,1,1/)
      LPIXELS=NAXES
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-flat=').GT.0)THEN
          FLATFILE=ARG(INDEX(ARG,'-flat=')+6:)
        ELSE IF(INDEX(ARG,'-dark=').GT.0)THEN
          DARKFILE=ARG(INDEX(ARG,'-dark=')+6:)
        ELSE IF(INDEX(ARG,'-fpixels=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-fpixels=')+9:),*)
     &      FPIXELS(1),FPIXELS(2),FPIXELS(3)
        ELSE IF(INDEX(ARG,'-lpixels=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-lpixels=')+9:),*)
     &      LPIXELS(1),LPIXELS(2),LPIXELS(3)
        ELSE IF(INDEX(ARG,'-o=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-o=')+3:),*) CRCTFILE
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          STOP
        END IF
      END DO
      WRITE(*,'(A)') ' target file: '//TRIM(TARGFILE)
      WRITE(*,'(A)') ' flat field file: '//TRIM(FLATFILE)
      WRITE(*,'(A)') ' dark current file: '//TRIM(DARKFILE)
      WRITE(*,'(A,I3,A,I3)') ' subimage size (width x height): ',
     &  LPIXELS(1)-FPIXELS(1)+1,' x',LPIXELS(2)-FPIXELS(2)+1
      WRITE(*,'(A,I5,A,I5)') ' from ',FPIXELS(3),' to ',LPIXELS(3)
      CALL FLATFIELDCORRECT(TARGFILE,FLATFILE,DARKFILE,
     &  FPIXELS,LPIXELS,CRCTFILE)
      STOP
      END PROGRAM FFC

      PROGRAM BMPTOFITS
C  Convert BMP picture to FITS file.
C  =================================
C
C  Usage:
C  ======
C  ./bmptofits bmpfile.bmp [fitsfile.fits]
C
      IMPLICIT NONE
      INTEGER :: STATUS,NARGS,BMPSIZ(2)
      CHARACTER(LEN=256) :: ARG,INFILE,OUTFILE,BASENAME,EXTNAME,PREFIX
      INTERFACE
      SUBROUTINE BMPSIZE(FILENAME,SIZ)
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: FILENAME
      INTEGER,INTENT(OUT) :: SIZ(2)
      END SUBROUTINE BMPSIZE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      END INTERFACE
      STATUS=0
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'kisa filename [-len-branch=n] [-num-branch=n]'//
     &    ' -psf=filename [-init=filename] [-prefix=output] [-n=numit]'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,INFILE)
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
      OUTFILE=TRIM(BASENAME)//'.fits'
      IF(NARGS .GT. 1)THEN
        CALL GET_COMMAND_ARGUMENT(2,OUTFILE)
      END IF
      CALL BMPSIZE(INFILE,BMPSIZ)
      STOP
      END PROGRAM BMPTOFITS
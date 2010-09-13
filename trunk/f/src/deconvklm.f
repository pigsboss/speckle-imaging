      PROGRAM MAIN
C  Kalman filter deconvolution executable.
C
C  Usage:
C  ======
C  deconvklm filename_g -psf=filename_h [-o=filename_f]
C
C  Purpose:
C  ========
C  g=conv(f,h). estimates f given g and h.
C
C  Arguments:
C  ==========
C  filename_g - input file. contains matrix g.
C  filename_h - input file. contains matrix h.
C  filename_f - output file. contains the estimate of matrix f.
C
C  Declaration:
C  ============
      IMPLICIT NONE
      INTEGER :: NAXES(3),NARGS,STATUS,K
      DOUBLE PRECISION, ALLOCATABLE :: DF(:,:),DG(:,:,:),DH(:,:)
      CHARACTER(LEN=256) :: ARG,FILEG,FILEH,FILEF,BASENAME,EXTNAME
C
      INTERFACE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE DECONVKLM(NX,NY,NFRMS,DG,DH,DF)
      INTEGER, INTENT(IN) :: NX,NY,NFRMS
      DOUBLE PRECISION, INTENT(IN) :: DG(NX,NY,NFRMS),DH(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DF(NX,NY)
      END SUBROUTINE DECONVKLM
      END INTERFACE
C
      STATUS=0
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,FILEG)
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'deconvklm filename_g -psf=filename -o=filename_f'
        STOP
      END IF
      CALL RESOLVEPATH(FILEG,BASENAME,EXTNAME)
      FILEF=TRIM(BASENAME)//'_klm.fits'
      FILEH=''
      CALL IMAGESIZE(FILEG,NAXES)
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-psf=').GT.0)THEN
          FILEH=ARG(INDEX(ARG,'-psf=')+5:)
        ELSE IF(INDEX(ARG,'-o=').GT.0)THEN
          FILEF=ARG(INDEX(ARG,'-o=')+3:)
        ELSE
          PRINT *,'unknown argument '//ARG
          STOP
        END IF
      END DO
      IF(LEN_TRIM(FILEH) .LT. 1)THEN
        PRINT *,'the PSF must be specified.'
        STOP
      END IF
      ALLOCATE(DG(NAXES(1),NAXES(2),NAXES(3)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        PRINT *,'out of memory.'
        STOP
      END IF
      ALLOCATE(DH(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        PRINT *,'out of memory.'
        STOP
      END IF
      CALL READIMAGE(FILEG,(/1,1,1/),(/NAXES(1),NAXES(2),NAXES(3)/),DG)
      CALL READIMAGE(FILEH,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DH)
      ALLOCATE(DF(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        PRINT *,'out of memory.'
        STOP
      END IF
      CALL DECONVKLM(NAXES(1),NAXES(2),NAXES(3),DG,DH,DF)
      DEALLOCATE(DG)
      DEALLOCATE(DH)
      DEALLOCATE(DF)
      STOP
      END PROGRAM MAIN
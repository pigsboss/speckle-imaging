      PROGRAM SPETOIMG
      IMPLICIT NONE
      INTEGER :: STATUS,NAXES(3),NARGS
      DOUBLE PRECISION, ALLOCATABLE :: DRHO(:,:),DPHI(:,:),DIMG(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZSPE(:,:)
      CHARACTER(LEN=256) :: MODFILE,PHAFILE,IMGFILE
C
      INTERFACE
      SUBROUTINE DFFTSHIFT(NX,NY,DX)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(INOUT) :: DX(NX,NY)
      END SUBROUTINE DFFTSHIFT
      SUBROUTINE DIFFTSHIFT(NX,NY,DX)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(INOUT) :: DX(NX,NY)
      END SUBROUTINE DIFFTSHIFT
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
      SUBROUTINE SPECTRUMTOIMAGE(NX,NY,ZSPE,DIMG)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(NX,NY)
      DOUBLE COMPLEX, INTENT(IN) :: ZSPE(NX,NY)
      END SUBROUTINE SPECTRUMTOIMAGE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      END INTERFACE
C
      STATUS=0
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS .LT. 3)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'spe2img modulus phase image'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,MODFILE)
      CALL GET_COMMAND_ARGUMENT(2,PHAFILE)
      CALL GET_COMMAND_ARGUMENT(3,IMGFILE)
C
      ALLOCATE(DRHO(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        STOP
      END IF
      ALLOCATE(DPHI(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        STOP
      END IF
      ALLOCATE(ZSPE(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        STOP
      END IF
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS .NE. 0)THEN
        WRITE(*,*)'error: out of memory.'
        STOP
      END IF
C
      CALL READIMAGE(MODFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DRHO)
      PRINT *,'input modulus: '//TRIM(MODFILE)
      CALL READIMAGE(PHAFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DPHI)
      PRINT *,'input phase: '//TRIM(PHAFILE)
      ZSPE=DCMPLX(DCOS(DPHI),DSIN(DPHI))*DRHO
      CALL SPECTRUMTOIMAGE(NAXES(1),NAXES(2),ZSPE,DIMG)
      CALL WRITEIMAGE(IMGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      PRINT *,'output image: '//TRIM(IMGFILE)
      DEALLOCATE(DRHO)
      DEALLOCATE(DPHI)
      DEALLOCATE(DIMG)
      DEALLOCATE(ZSPE)
      STOP
      END PROGRAM SPETOIMG
      PROGRAM BGFIT
C  Usage:
C  ======
C  bgfit filename_avg filename_bg [-r=radius] [-fit=fit_method]
C
C  Purpose:
C  ========
C  Estimate background parameters and background fitting the edge of the given
C  average image.
C
C  Arguments:
C  ==========
C  filename_avg - filename of average image.
C  filename_bg  - filename of background, output.
C  radius       - radius of the signal region.
C  fit_method   - fitting method, such as 'p0', 'p2', and 'p4'.
C
      IMPLICIT NONE
      INTEGER :: NAXES(3),NPARAMS,K,NARGS
      DOUBLE PRECISION :: DR
      DOUBLE PRECISION, ALLOCATABLE :: DIMG(:,:),DBG(:,:),DPARAMS(:)
      CHARACTER(LEN=256) :: AVGFILE,BGFILE,FTMETHOD,ARG
      INTERFACE
      SUBROUTINE BGFIT2P0(M,N,DR,DIMG,DBG,DPARAMS)
      INTEGER, INTENT(IN) :: M,N
      DOUBLE PRECISION, INTENT(IN) :: DR,DIMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(M,N),DPARAMS(*)
      END SUBROUTINE BGFIT2P0
      SUBROUTINE BGFIT2P2(M,N,DR,DIMG,DBG,DPARAMS)
      INTEGER, INTENT(IN) :: M,N
      DOUBLE PRECISION, INTENT(IN) :: DR,DIMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(M,N),DPARAMS(*)
      END SUBROUTINE BGFIT2P2
      SUBROUTINE BGFIT2P4(M,N,DR,DIMG,DBG,DPARAMS)
      INTEGER, INTENT(IN) :: M,N
      DOUBLE PRECISION, INTENT(IN) :: DR,DIMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(M,N),DPARAMS(*)
      END SUBROUTINE BGFIT2P4
      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER*(*) :: IMGFILE
      END SUBROUTINE IMAGESIZE
      END INTERFACE
      CALL GET_COMMAND_ARGUMENT(1,AVGFILE)
      CALL GET_COMMAND_ARGUMENT(2,BGFILE)
      CALL IMAGESIZE(AVGFILE,NAXES)
      DR=0.5*DBLE(MIN(NAXES(1),NAXES(2)))
      FTMETHOD='p4'
      DO K=3,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-r=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-r=')+3:256),*) DR
        ELSE IF(INDEX(ARG,'-fit=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-fit=')+5:256),*) FTMETHOD
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          RETURN
        END IF
      END DO
      ALLOCATE(DIMG(NAXES(2),NAXES(1)))
      ALLOCATE(DBG(NAXES(2),NAXES(1)))
      ALLOCATE(DPARAMS(NAXES(2)*NAXES(1)))
      CALL READIMAGE(AVGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      IF(FTMETHOD .EQ. 'p0')THEN
        NPARAMS=1
        CALL BGFIT2P0(NAXES(2),NAXES(1),DR,DIMG,DBG,DPARAMS)
      ELSE IF(FTMETHOD .EQ. 'p2')THEN
        NPARAMS=6
        CALL BGFIT2P2(NAXES(2),NAXES(1),DR,DIMG,DBG,DPARAMS)
      ELSE IF(FTMETHOD .EQ. 'p4')THEN
        NPARAMS=15
        CALL BGFIT2P4(NAXES(2),NAXES(1),DR,DIMG,DBG,DPARAMS)
      ELSE
        PRINT *,'Unknown fitting method '//TRIM(FTMETHOD)
        RETURN
      END IF
      CALL WRITEIMAGE(BGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
      DO K=1,NPARAMS
        WRITE(ARG,*)K
        WRITE(*,'(A,ES10.3)')' a_'//TRIM(ADJUSTL(ARG))//' = ',
     &    DPARAMS(K)
      END DO
      DEALLOCATE(DIMG)
      DEALLOCATE(DBG)
      DEALLOCATE(DPARAMS)
      STOP
      END PROGRAM BGFIT

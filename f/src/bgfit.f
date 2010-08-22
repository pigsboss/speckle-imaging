      PROGRAM BGFIT
C  Usage:
C  ======
C  bgfit filename_avg [-r=radius] [-fit=fit_method] [-prefix=...]
C
C  Purpose:
C  ========
C  Estimate background parameters and background fitting the edge of the given
C  average image.
C
C  Arguments:
C  ==========
C  filename_avg - filename of average image.
C  radius       - radius of the signal region.
C  fit_method   - fitting method, such as 'p0', 'p2', and 'p4'.
C  prefix       - prefix of output filenames.
C
      IMPLICIT NONE
      INTEGER :: NAXES(3),NPARAMS,K,X,Y,NARGS,N,MAXN
      DOUBLE PRECISION :: DR,DXC,DYC,DTMP
      DOUBLE PRECISION, ALLOCATABLE :: DIMG(:,:),DBG(:,:),DFLAG(:,:),
     &  DPARAMS(:)
      CHARACTER(LEN=256) :: AVGFILE,PREFIX,FTMETHOD,ARG,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER*(*), INTENT(IN) :: PATH
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
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
      SUBROUTINE BGFIT2PN(NX,NY,DR,DIMG,N,DBG,DB)
      INTEGER, INTENT(IN) :: NX,NY,N
      DOUBLE PRECISION, INTENT(IN) :: DR,DIMG(NY,NX)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(NY,NX),DB(NY*NX)
      END SUBROUTINE BGFIT2PN
      END INTERFACE
      MAXN=8
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,AVGFILE)
      CALL IMAGESIZE(AVGFILE,NAXES)
      DR=0.5*DBLE(MIN(NAXES(1),NAXES(2)))
      FTMETHOD='all'
      CALL RESOLVEPATH(AVGFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-r=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-r=')+3:),*) DR
        ELSE IF(INDEX(ARG,'-fit=').GT.0)THEN
          FTMETHOD=ARG(INDEX(ARG,'-fit=')+5:)
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=ARG(INDEX(ARG,'-prefix=')+8:)
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          RETURN
        END IF
      END DO
      ALLOCATE(DIMG(NAXES(2),NAXES(1)))
      ALLOCATE(DFLAG(NAXES(2),NAXES(1)))
      ALLOCATE(DBG(NAXES(2),NAXES(1)))
      ALLOCATE(DPARAMS(NAXES(2)*NAXES(1)))
      DXC=0.5D0*DBLE(NAXES(1)+1)
      DYC=0.5D0*DBLE(NAXES(2)+1)
      DO X=1,NAXES(1)
        DO Y=1,NAXES(2)
          DTMP=DSQRT((DBLE(X)-DXC)*(DBLE(X)-DXC)+
     &      (DBLE(Y)-DYC)*(DBLE(Y)-DYC))
          IF(DTMP.GE.DR)THEN
            DFLAG(Y,X)=1.0D0
          ELSE
            DFLAG(Y,X)=0.0D0
          END IF
        END DO
      END DO
      CALL READIMAGE(AVGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      DBG=DIMG*DFLAG
      CALL WRITEIMAGE(TRIM(PREFIX)//'_bg.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'defined background: ',TRIM(PREFIX)//'_bg.fits'
      IF(FTMETHOD .EQ. 'all')THEN
        DO N=0,MAXN
          WRITE(FTMETHOD,*)N
          FTMETHOD=TRIM(ADJUSTL(FTMETHOD))
          PRINT *,'using '//TRIM(FTMETHOD)//'-th polynomials:'
          FTMETHOD='p'//TRIM(ADJUSTL(FTMETHOD))
          CALL BGFIT2PN(NAXES(1),NAXES(2),DR,DIMG,N,DBG,DPARAMS)
          NPARAMS=(N+1)*(N+2)/2
          DO K=1,NPARAMS
            WRITE(ARG,*)K
            WRITE(*,'(A,ES10.3)')' a_'//TRIM(ADJUSTL(ARG))//' = ',
     &        DPARAMS(K)
          END DO
          CALL WRITEIMAGE(TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_bg.fits'
     &      ,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
          PRINT *,'fitted background: ',
     &      TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_bg.fits'
          DBG=(DIMG-DBG)*DFLAG
          CALL WRITEIMAGE(TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_res.fits'
     &      ,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
          PRINT *,'residual: ',
     &      TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_res.fits'
        END DO
      ELSE
        READ(FTMETHOD(2:),*) N
        CALL BGFIT2PN(NAXES(1),NAXES(2),DR,DIMG,N,DBG,DPARAMS)
        NPARAMS=(N+1)*(N+2)/2
        DO K=1,NPARAMS
          WRITE(ARG,*)K
          WRITE(*,'(A,ES10.3)')' a_'//TRIM(ADJUSTL(ARG))//' = ',
     &      DPARAMS(K)
        END DO
        CALL WRITEIMAGE(TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_bg.fits',
     &    (/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
        PRINT *,'fitted background: ',
     &    TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_bg.fits'
        DBG=(DIMG-DBG)*DFLAG
        CALL WRITEIMAGE(TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_res.fits',
     &    (/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
        PRINT *,'residual: ',
     &    TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_res.fits'
      END IF
      DEALLOCATE(DIMG)
      DEALLOCATE(DBG)
      DEALLOCATE(DPARAMS)
      DEALLOCATE(DFLAG)
      STOP
      END PROGRAM BGFIT

      PROGRAM BGFIT
C  Usage:
C  ======
C  bgfit filename_avg [-flag=filename] [-round=x_c,y_c,radius]
C    [-fit=fit_method] [-prefix=output]
C
C  Purpose:
C  ========
C  Estimate background parameters and background fitting the edge of the given
C  average image.
C
C  Arguments:
C  ==========
C
      IMPLICIT NONE
      INTEGER :: STATUS,NAXES(3),NPARAMS,K,X,Y,NARGS,N,MAXN,NSAMPLES
      DOUBLE PRECISION :: DR,DXC,DYC,DTMP
      DOUBLE PRECISION, ALLOCATABLE :: DIMG(:,:),DBG(:,:),DFLAG(:,:),
     &  DPARAMS(:)
      CHARACTER(LEN=256) :: AVGFILE,PREFIX,FTMETHOD,ARG,BASENAME,
     &  EXTNAME,FLAGFILE
      INTERFACE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER*(*), INTENT(IN) :: PATH
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER*(*) :: IMGFILE
      END SUBROUTINE IMAGESIZE
      SUBROUTINE BGFIT2PN(NX,NY,DFLAG,DIMG,N,DBG,DB)
      INTEGER, INTENT(IN) :: NX,NY,N
      DOUBLE PRECISION, INTENT(IN) :: DFLAG(NX,NY),DIMG(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(NX,NY),DB(*)
      END SUBROUTINE BGFIT2PN
      END INTERFACE
      STATUS=0
      MAXN=8
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'bgfit filename_avg [-flag=filename]'//
     &    ' [-round=x_c,y_c,radius] [-fit=fit_method]'//
     &    ' [-prefix=output]'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,AVGFILE)
      CALL IMAGESIZE(AVGFILE,NAXES)
      DR=0.5*DBLE(MIN(NAXES(1),NAXES(2)))
      DXC=0.5D0*DBLE(NAXES(1)+1)
      DYC=0.5D0*DBLE(NAXES(2)+1)
      FTMETHOD='all'
      FLAGFILE=''
      CALL RESOLVEPATH(AVGFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)
      ALLOCATE(DFLAG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      DO X=1,NAXES(1)
        DO Y=1,NAXES(2)
          DTMP=DSQRT((DBLE(X)-DXC)*(DBLE(X)-DXC)+
     &      (DBLE(Y)-DYC)*(DBLE(Y)-DYC))
          IF(DTMP.GE.DR)THEN
            DFLAG(X,Y)=1.0D0
          ELSE
            DFLAG(X,Y)=0.0D0
          END IF
        END DO
      END DO
      NSAMPLES=NINT(SUM(DFLAG))
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-round=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-round=')+7:),*)DXC,DYC,DR
          DO X=1,NAXES(1)
            DO Y=1,NAXES(2)
              DTMP=DSQRT((DBLE(X)-DXC)*(DBLE(X)-DXC)+
     &          (DBLE(Y)-DYC)*(DBLE(Y)-DYC))
              IF(DTMP.GE.DR)THEN
                DFLAG(X,Y)=1.0D0
              ELSE
                DFLAG(X,Y)=0.0D0
              END IF
            END DO
          END DO
          NSAMPLES=NINT(SUM(DFLAG))
        ELSE IF(INDEX(ARG,'-flag=').GT.0)THEN
          FLAGFILE=ARG(INDEX(ARG,'-flag=')+6:)
          CALL READIMAGE(FLAGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),
     &      DFLAG)
        ELSE IF(INDEX(ARG,'-fit=').GT.0)THEN
          FTMETHOD=ARG(INDEX(ARG,'-fit=')+5:)
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=ARG(INDEX(ARG,'-prefix=')+8:)
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          STOP
        END IF
      END DO
      CALL WRITEIMAGE(TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_flag.fits'
     &  ,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DFLAG)
      PRINT *,'defined background region: ',
     &  TRIM(PREFIX)//'_'//TRIM(FTMETHOD)//'_flag.fits'
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(AVGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      ALLOCATE(DBG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      DBG=DIMG*DFLAG
      CALL WRITEIMAGE(TRIM(PREFIX)//'_bg.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'defined background: ',TRIM(PREFIX)//'_bg.fits'
      ALLOCATE(DPARAMS(NAXES(1)*NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      IF(FTMETHOD .EQ. 'all')THEN
        DO N=0,MAXN
          WRITE(FTMETHOD,*)N
          FTMETHOD=TRIM(ADJUSTL(FTMETHOD))
          PRINT *,'using '//TRIM(FTMETHOD)//'-th polynomials:'
          FTMETHOD='p'//TRIM(ADJUSTL(FTMETHOD))
          CALL BGFIT2PN(NAXES(1),NAXES(2),DFLAG,DIMG,N,DBG,DPARAMS)
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
        CALL BGFIT2PN(NAXES(1),NAXES(2),DFLAG,DIMG,N,DBG,DPARAMS)
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

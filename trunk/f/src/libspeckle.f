      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      IMPLICIT NONE
      INTEGER :: STATUS,UNIT,RWMODE,BLOCKSIZE,NAXIS,NFOUND,GROUP
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      RWMODE=0
      CALL FTOPEN(UNIT,IMGFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'IMAGESIZE failed to read the NAXIS keywords.'
        RETURN
      ENDIF
      CALL FTCLOS(UNIT, STATUS)
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE IMAGESIZE
C ******************************************************************************
      SUBROUTINE SUBIMAGE(IMGFILE,FPIXELS,LPIXELS,SUBFILE,BUFFERSIZE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      CHARACTER*(*), INTENT(IN) :: IMGFILE,SUBFILE
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      INTEGER :: STATUS,K,M,N,NPIXELS,NFRAMES,LBUFFER,NBUFFER
      INTEGER :: L,L1,L2
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:)
      STATUS=0
      CALL DELETEFILE(SUBFILE,STATUS)
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=M*N
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      IF (PRESENT(BUFFERSIZE))THEN
        LBUFFER=INT(FLOOR(DBLE(BUFFERSIZE*1024*1024/8)/DBLE(NPIXELS)))
      ELSE
        LBUFFER=INT(FLOOR(DBLE(32*1024*1024/8)/DBLE(NPIXELS)))
      END IF
      NBUFFER=INT(CEILING(DBLE(NFRAMES)/DBLE(LBUFFER)))
      ALLOCATE(DBUF(M,N,LBUFFER))
      DO K=1,NBUFFER
        L1=(K-1)*LBUFFER+1
        L2=MIN(K*LBUFFER,NFRAMES)
C       PRINT *,'READIMAGE from',L1,'to',L2
        CALL READIMAGE(IMGFILE,(/FPIXELS(1),FPIXELS(2),L1/),
     &    (/LPIXELS(1),LPIXELS(2),L2/),DBUF)
        CALL APPENDIMAGE(SUBFILE,(/1,1,L1/),(/M,N,L2/),DBUF)
      END DO
      DEALLOCATE(DBUF)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE SUBIMAGE
C ******************************************************************************
      SUBROUTINE FLATFIELDCORRECT(TARGFILE,FLATFILE,DARKFILE,
     &  FPIXELS,LPIXELS,CRCTFILE)
C  Flat-field correction subroutine.
C
C  Usage:
C  =====
C  Input direct observation results, flat field file, as well as dark current
C  file. Return flat-field corrected observation results as output.
C
C                 I_original - I_dark
C  I_corrected = ---------------------
C                   I_flat - I_dark
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: STATUS,UNIT,RWMODE,BLOCKSIZE,GROUP
      INTEGER :: K,M,N,NPIXELS,NFRAMES,NAXES(3),NFOUND
      DOUBLE PRECISION, ALLOCATABLE :: DFLAT(:),DDARK(:),DTARG(:),
     &  DCRCT(:)
      CHARACTER*(*), INTENT(IN) :: TARGFILE,FLATFILE,DARKFILE,CRCTFILE
      CHARACTER(LEN=256) :: KEYVAL,COMMENT
      INTERFACE
      SUBROUTINE AVERAGE(AVGFILE,FPIXELS,LPIXELS,DAVG,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER, INTENT(INOUT), OPTIONAL :: BUFFERSIZE
      DOUBLE PRECISION :: DAVG(LPIXELS(1)+1-FPIXELS(1),
     &  LPIXELS(2)+1-FPIXELS(2))
      CHARACTER*(*), INTENT(IN) :: AVGFILE
      END SUBROUTINE AVERAGE
      END INTERFACE
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      RWMODE=0
      CALL FTOPEN(UNIT,TARGFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(TARGFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      M=NAXES(1)
      N=NAXES(2)
      CALL FTCLOS(UNIT, STATUS)
      CALL FTOPEN(UNIT,FLATFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(FLATFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      IF ((NAXES(1).NE.M) .OR. (NAXES(2).NE.N))THEN
        PRINT *,TRIM(TARGFILE)//' and '//TRIM(FLATFILE)//
     &    ' differ in size.'
        RETURN
      END IF
      CALL FTCLOS(UNIT, STATUS)
      ALLOCATE(DFLAT(NAXES(1)*NAXES(2)))
      CALL AVERAGE(FLATFILE,(/FPIXELS(1),FPIXELS(2),1/),
     &  (/LPIXELS(1),LPIXELS(2),NAXES(3)/),DFLAT)
      CALL FTOPEN(UNIT,DARKFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(DARKFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      IF ((NAXES(1).NE.M) .OR. (NAXES(2).NE.N))THEN
        PRINT *,TRIM(TARGFILE)//' and '//TRIM(DARKFILE)//
     &    ' differ in size.'
        RETURN
      END IF
      CALL FTCLOS(UNIT, STATUS)
      ALLOCATE(DDARK(NAXES(1)*NAXES(2)))
      CALL AVERAGE(DARKFILE,(/FPIXELS(1),FPIXELS(2),1/),
     &  (/LPIXELS(1),LPIXELS(2),NAXES(3)/),DDARK)
      IF (MAXVAL(DDARK) .GE. MINVAL(DFLAT))THEN
        PRINT *,'Maximum of dark current larger than or equal'
     &    //' to minimum of flat field.'
        RETURN
      END IF
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=M*N
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      ALLOCATE(DTARG(NPIXELS))
      ALLOCATE(DCRCT(NPIXELS))
      CALL DELETEFILE(CRCTFILE,STATUS)
      DO K=FPIXELS(3),LPIXELS(3)
        CALL READIMAGE(TARGFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &    (/LPIXELS(1),LPIXELS(2),K/),DTARG)
        DCRCT=(DTARG-DDARK)/(DFLAT-DDARK)
        CALL APPENDIMAGE(CRCTFILE,(/1,1,K/),(/M,N,K/),DCRCT)
      END DO
      DEALLOCATE(DTARG)
      DEALLOCATE(DFLAT)
      DEALLOCATE(DDARK)
      DEALLOCATE(DCRCT)
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE FLATFIELDCORRECT
C ******************************************************************************
      SUBROUTINE PRINTERROR(STATUS)
      INTEGER, INTENT(IN) :: STATUS
      CHARACTER :: ERRTEXT*30,ERRMESSAGE*80
      IF (STATUS .LE. 0)RETURN
      CALL FTGERR(STATUS,ERRTEXT)
      PRINT *,'FITSIO ERROR STATUS =',STATUS,': ',ERRTEXT
      CALL FTGMSG(ERRMESSAGE)
      DO WHILE (ERRMESSAGE .NE. ' ')
        PRINT *,ERRMESSAGE
        CALL FTGMSG(ERRMESSAGE)
      ENDDO
      RETURN
      END SUBROUTINE PRINTERROR
C ******************************************************************************
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      INTEGER, INTENT(INOUT) :: STATUS
      INTEGER :: UNIT,BLOCKSIZE
      CHARACTER*(*), INTENT(IN) :: FILENAME
      IF (STATUS .GT. 0)RETURN
      CALL FTGIOU(UNIT,STATUS)
      CALL FTOPEN(UNIT,FILENAME,1,BLOCKSIZE,STATUS)
      IF (STATUS .EQ. 0)THEN
        CALL FTDELT(UNIT,STATUS)
      ELSE IF (STATUS .EQ. 103)THEN
        STATUS=0
        CALL FTCMSG
      ELSE
        STATUS=0
        CALL FTCMSG
        CALL FTDELT(UNIT,STATUS)
      END IF
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE DELETEFILE
C ******************************************************************************
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER :: STATUS,UNIT,BLOCKSIZE,BITPIX,NAXIS,GROUP
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: NAXES(3)
      DOUBLE PRECISION, INTENT(IN) :: ARRAY(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      BLOCKSIZE=1
      BITPIX=-64
      NAXIS=3
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      CALL DELETEFILE(FILENAME,STATUS)
      CALL FTINIT(UNIT,FILENAME,BLOCKSIZE,STATUS)
      CALL FTPHPS(UNIT,BITPIX,NAXIS,NAXES,STATUS)
      GROUP=1
      CALL FTPSSD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,ARRAY,STATUS)
      CALL FTCLOS(UNIT,STATUS)
      CALL FTFIOU(UNIT,STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE WRITEIMAGE
C ******************************************************************************
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: STATUS,UNIT,READWRITE,BLOCKSIZE,NAXIS,NFOUND
      INTEGER :: GROUP,INCS(3),NAXES(3)
      DOUBLE PRECISION, INTENT(OUT) :: ARRAY(*)
      DOUBLE PRECISION :: NULLVAL
      LOGICAL :: ANYF
      CHARACTER*(*), INTENT(IN) :: FILENAME
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      READWRITE=0
      CALL FTOPEN(UNIT,FILENAME,READWRITE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'READIMAGE failed to read the NAXIS keywords.'
        RETURN
      ENDIF
      GROUP=1
      NAXIS=3
      INCS(1)=1
      INCS(2)=1
      INCS(3)=1
      NULLVAL=-999
      CALL FTGSVD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,INCS,
     &  NULLVAL,ARRAY,ANYF,STATUS)
      IF (ANYF)THEN
        PRINT *,'One or more pixels are undefined in the FITS image.'
        RETURN
      ENDIF
      CALL FTCLOS(UNIT, STATUS)
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE READIMAGE
C ******************************************************************************
      SUBROUTINE AVERAGE(FILENAME,FPIXELS,LPIXELS,ARRAY,BUFFERSIZE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER, INTENT(INOUT), OPTIONAL :: BUFFERSIZE
      INTEGER :: L,LF,LT,LR,K,NPIXELS,M,N,NFRAMES,LBUFFER
      CHARACTER*(*), INTENT(IN) :: FILENAME
      DOUBLE PRECISION, INTENT(OUT) :: ARRAY(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(2)-FPIXELS(2)+1)
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:,:,:)
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=M*N
      IF (PRESENT(BUFFERSIZE))THEN
        LBUFFER=INT(FLOOR(DBLE(BUFFERSIZE)*1024*128/DBLE(NPIXELS)))
      ELSE
        LBUFFER=INT(FLOOR(DBLE(4*1024*1024)/DBLE(NPIXELS)))
      END IF
      ALLOCATE(BUFFER(M,N,LBUFFER))
      ARRAY=0
      LF=1
      DO L=1,INT(CEILING(DBLE(NFRAMES)/DBLE(LBUFFER)))
        LT=MIN(LPIXELS(3),LF+LBUFFER-1)
        CALL READIMAGE(FILENAME,(/FPIXELS(1),FPIXELS(2),LF/),
     &    (/LPIXELS(1),LPIXELS(2),LT/),BUFFER(1:M,1:N,1:LBUFFER))
        LR=LT-LF+1
        LF=LT+1
        DO K=1,LR
          ARRAY=ARRAY+BUFFER(1:M,1:N,K)
        END DO
      END DO
      ARRAY=ARRAY/DBLE(NFRAMES)
      DEALLOCATE(BUFFER)
      RETURN
      END SUBROUTINE AVERAGE
C ******************************************************************************
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER*(*), INTENT(IN) :: PATH
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      CHARACTER :: FILENAME*80
      INTEGER :: K,L
      K=SCAN(PATH,'/',.TRUE.)
      L=LEN_TRIM(PATH)
      READ(PATH(K+1:L),*) FILENAME
      K=SCAN(FILENAME,'.')
      L=LEN_TRIM(FILENAME)
      READ(FILENAME(1:K-1),*) BASENAME
      K=SCAN(FILENAME,'.',.TRUE.)
      READ(FILENAME(K+1:L),*) EXTNAME
      RETURN
      END SUBROUTINE RESOLVEPATH
C ******************************************************************************
      SUBROUTINE BGFIT2P0(M,N,D,IMG,BG,B)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      NPARAMS=1
      NSAMPLES=0
      B(1)=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
            B(1)=B(1)+IMG(I,J)
          END IF
        END DO
      END DO
      B(1)=B(1)/DBLE(NSAMPLES)
      BG=B(1)
      RETURN
      END SUBROUTINE BGFIT2P0
C ******************************************************************************
      SUBROUTINE BGFIT2P2(M,N,D,IMG,BG,B)
C  2-dimensional background fitting subroutine using 2nd polynomials.
C
C  Purpose
C  =======
C  z = a_0 + a_1*x + a_2*y + a_3*x^2 + a_4*x*y + a_5*y^2
C  x is column number.
C  y is row number.
C  z is value of IMG at (x,y).
C  Try to determine parameters a_i, i=0,1,2,3,4,5.
C
C  Arguments
C  =========
C  M is the number of rows of matrix IMG.
C  N is the number of columns of matrix IMG.
C  IMG is the matrix the background of which is to be fitted.
C  BG,B are the output of this subroutine.
C  Each row of A is vector (1, x, y, x^2, x*y, y^2) for a specific (x,y).
C  B is (z_1, z_2, z_3, ..., z_n)'. The subscript of z denotes different
C  location.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:)
      NPARAMS=6
      NSAMPLES=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(A(NSAMPLES,6))
      K=1
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            A(K,:)=(/DBLE(1),X,Y,X*X,X*Y,Y*Y/)
            B(K)=IMG(INT(Y),INT(X))
            K=K+1
          ENDIF
        ENDDO
      ENDDO
      LDA=NSAMPLES
      LDB=MAX(NSAMPLES,6)
      LWORK=-1
      ALLOCATE(WORK(1))
      CALL DGELS('N',NSAMPLES,6,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      LWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DGELS('N',NSAMPLES,6,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      IF (INFO.GT.0)THEN
        PRINT *,'The least-squares solution could not be computed ',
     &    'because A does not have full rank.'
        RETURN
      ENDIF
      IF (INFO.LT.0)THEN
        PRINT *,'The argument has illegal value: ',ABS(INFO)
        RETURN
      ENDIF
      DO I=1,M
        DO J=1,N
          BG(I,J)=B(1)+B(2)*DBLE(J)+B(3)*DBLE(I)+B(4)*DBLE(J)*DBLE(J)
     &      +B(5)*DBLE(J)*DBLE(I)+B(6)*DBLE(I)*DBLE(I)
        ENDDO
      ENDDO
      DEALLOCATE(WORK)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE BGFIT2P2
C ******************************************************************************
      SUBROUTINE BGFIT2P4(M,N,D,IMG,BG,B)
C  2-dimensional background fitting subroutine using 4th polynomials.
C
C  Purpose
C  =======
C  z = a_0 + 
C    a_1*x + a_2*y + 
C    a_3*x^2 + a_4*x*y + a_5*y^2 + 
C    a_6*x^3 + a_7*x^2*y + a_8*x*y^2 + a_9*y^3 +
C    a_10*x^4 + a_11*x^3*y + a_12*x^2*y^2 + a_13*x*y^3 + a_14*y^4
C  x is column number.
C  y is row number.
C  z is value of IMG at (x,y).
C  Try to determine parameters a_i, i=0,1,2,3,...,14.
C
C  Arguments
C  =========
C  M is the number of rows of matrix IMG.
C  N is the number of columns of matrix IMG.
C  D is distance. Elements whose distance from the centre of IMG is larger than
C    D are sampled. BG is fitted to the sampled elements.
C  IMG is the matrix the background of which is to be fitted.
C  BG is the output of this subroutine.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:)
      NPARAMS=15
      NSAMPLES=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(A(NSAMPLES,NPARAMS))
      K=1
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (DSQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            A(K,:)=(/DBLE(1),X,Y,X*X,X*Y,Y*Y,
     &        X*X*X,X*X*Y,X*Y*Y,Y*Y*Y,
     &        X*X*X*X,X*X*X*Y,X*X*Y*Y,X*Y*Y*Y,Y*Y*Y*Y/)
            B(K)=IMG(INT(Y),INT(X))
            K=K+1
          ENDIF
        ENDDO
      ENDDO
      LDA=NSAMPLES
      LDB=MAX(NSAMPLES,NPARAMS)
      LWORK=-1
      ALLOCATE(WORK(1))
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      LWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      IF (INFO.GT.0)THEN
        PRINT *,'The least-squares solution could not be computed ',
     &    'because A does not have full rank.'
        RETURN
      ENDIF
      IF (INFO.LT.0)THEN
        PRINT *,'The argument has illegal value: ',ABS(INFO)
        RETURN
      ENDIF
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          BG(I,J)=B(1)+B(2)*X+B(3)*Y+
     &      B(4)*X*X+B(5)*X*Y+B(6)*Y*Y+
     &      B(7)*X*X*X+B(8)*X*X*Y+B(9)*X*Y*Y+B(10)*Y*Y*Y+
     &      B(11)*X*X*X*X+B(12)*X*X*X*Y+B(13)*X*X*Y*Y+B(14)*X*Y*Y*Y+
     &      B(15)*Y*Y*Y*Y
        ENDDO
      ENDDO
      DEALLOCATE(WORK)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE BGFIT2P4
C ******************************************************************************
      SUBROUTINE GETSNR(M,N,D,SNR,IMG,FTMETHOD)
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN1D,PLAN2D
      INTEGER :: M,N,I,J,K,NPIXS,NSPLS,INFO
      DOUBLE PRECISION :: SNR,D,X,Y,XC,YC,PS,PN
      DOUBLE PRECISION :: IMG(M,N),CLN(M,N),WORK(M,N)
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:)
      DOUBLE COMPLEX :: ZIN2D(M,N),ZOUT2D(M,N)
      DOUBLE COMPLEX, ALLOCATABLE :: ZIN1D(:),ZOUT1D(:)
      CHARACTER*(*) :: FTMETHOD
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
        RETURN
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      NPIXS=M*N
      NSPLS=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSPLS=NSPLS+1
          END IF
        END DO
      END DO
      ALLOCATE(BUFFER(NSPLS))
      ALLOCATE(ZIN1D(NSPLS))
      ALLOCATE(ZOUT1D(NSPLS))
      CALL DFFTW_PLAN_DFT_1D(PLAN1D,NSPLS,ZIN1D,ZOUT1D,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLAN2D,M,N,ZIN2D,ZOUT2D,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      IF (FTMETHOD .EQ. 'P2') THEN
        CALL BGFIT2P2(M,N,D,IMG,WORK,BUFFER)
      ELSE IF (FTMETHOD .EQ. 'P4') THEN
        CALL BGFIT2P4(M,N,D,IMG,WORK,BUFFER)
      ELSE IF (FTMETHOD .EQ. 'P0') THEN
        CALL BGFIT2P0(M,N,D,IMG,WORK,BUFFER(1))
      ELSE
        PRINT *,' Unknown fitting method.'
        RETURN
      END IF
      CLN=IMG-WORK
      K=0
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            K=K+1
            BUFFER(K)=CLN(I,J)
          END IF
        END DO
      END DO
      ZIN2D(1:M,1:N)=CMPLX(CLN(1:M,1:N))
      CALL DFFTW_EXECUTE_DFT(PLAN2D,ZIN2D,ZOUT2D)
      PS=DZNRM2(M*N,RESHAPE(ZOUT2D,(/NPIXS,1/)),1)
      PS=PS*PS/DBLE(NPIXS)/SQRT(DBLE(NPIXS))
      ZIN1D(1:NSPLS)=CMPLX(BUFFER(1:NSPLS))
      CALL DFFTW_EXECUTE_DFT(PLAN1D,ZIN1D,ZOUT1D)
      PN=DZNRM2(NSPLS,ZOUT1D(1:NSPLS),1)
      PN=PN*PN/DBLE(NSPLS)/SQRT(DBLE(NSPLS))
      SNR=PS/PN
      CALL DFFTW_DESTROY_PLAN(PLAN1D)
      CALL DFFTW_DESTROY_PLAN(PLAN2D)
C     CALL DFFTW_CLEANUP_THREADS()
      DEALLOCATE(ZIN1D)
      DEALLOCATE(BUFFER)
      DEALLOCATE(ZOUT1D)
      RETURN
      END SUBROUTINE GETSNR
C ******************************************************************************
      SUBROUTINE DECONVWNR(M,N,DG,DF,DH,DSNR)
C  Wiener deconvolution subroutine.
C
C  Purpose:
C  ========
C  G = conv(F, H). This routine returns F.
C
C  Arguments:
C  ==========
C  M - Number of rows of G.
C  N - Number of columns of G.
C  DG - Input.
C  DF - Output. size(F) = size(G).
C  DH - Input. size(PSF) = size (G)
C  DSNR - Signal-to-noise ratio.
C 
C  Declarations:
C  =============
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER :: M,N,INFO
      DOUBLE PRECISION :: DG(M,N),DF(M,N),DH(M,N)
      DOUBLE PRECISION :: DSNR
      DOUBLE COMPLEX :: ZG(M,N),ZH(M,N)
      DOUBLE COMPLEX :: ZIN(M,N),ZOUT(M,N),ZDECONV(M,N)
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=CMPLX(DG)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZG=ZOUT
      ZIN=CMPLX(DH)
      CALL ZFFTSHIFT(M,N,ZIN)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZH=ZOUT
      ZDECONV=CONJG(ZH)/(ZH*CONJG(ZH)+CMPLX(DBLE(1)/DSNR))
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZG*ZDECONV
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DF=DBLE(ZOUT)
      CALL DFFTW_DESTROY_PLAN(PLAN)
C     CALL DFFTW_CLEANUP_THREADS()
      RETURN
      END SUBROUTINE DECONVWNR
C ******************************************************************************
      SUBROUTINE ZFFTSHIFT(M,N,ZX)
C  Shift zero-frequency component to the centre of spectrum.
      INTEGER, INTENT(IN) :: M,N
      DOUBLE COMPLEX, INTENT(INOUT) :: ZX(M,N)
      ZX=CSHIFT(CSHIFT(ZX,INT(FLOOR(0.5*DBLE(M))),1),
     &  INT(FLOOR(0.5*DBLE(N))),2)
      RETURN
      END SUBROUTINE ZFFTSHIFT
C ******************************************************************************
      SUBROUTINE DFFTSHIFT(M,N,DX)
C  Shift zero-frequency component to the centre of spectrum.
      INTEGER, INTENT(IN) :: M,N
      DOUBLE PRECISION, INTENT(INOUT) :: DX(M,N)
      DX=CSHIFT(CSHIFT(DX,INT(FLOOR(0.5*DBLE(M))),1),
     &  INT(FLOOR(0.5*DBLE(N))),2)
      RETURN
      END SUBROUTINE DFFTSHIFT
C ******************************************************************************
      SUBROUTINE ZIFFTSHIFT(M,N,ZX)
C  Shift zero-frequency component to (1,1) position of spectrum.
      INTEGER, INTENT(IN) :: M,N
      DOUBLE COMPLEX, INTENT(INOUT) :: ZX(M,N)
      ZX=CSHIFT(CSHIFT(ZX,INT(FLOOR(-0.5*DBLE(M))),1),
     &  INT(FLOOR(-0.5*DBLE(N))),2)
      RETURN
      END SUBROUTINE ZIFFTSHIFT
C ******************************************************************************
      SUBROUTINE DIFFTSHIFT(M,N,DX)
C  Shift zero-frequency component to (1,1) position of spectrum.
      INTEGER, INTENT(IN) :: M,N
      DOUBLE PRECISION, INTENT(INOUT) :: DX(M,N)
      DX=CSHIFT(CSHIFT(DX,INT(FLOOR(-0.5*DBLE(M))),1),
     &  INT(FLOOR(-0.5*DBLE(N))),2)
      RETURN
      END SUBROUTINE DIFFTSHIFT
C ******************************************************************************
      SUBROUTINE DCORR2D(M,N,DA,DB,DC)
C  correlation coefficients of 2-dimensional array.
C
C  Purpose:
C  ========
C  c = corr(a, b). calculates c given a and b using fast fourier transform.
C  fft(c) = fft(a) * conj(fft(b)) (periodic extension is implied).
C
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER :: M,N
      DOUBLE PRECISION, DIMENSION(M,N) :: DA,DB,DC
      DOUBLE COMPLEX, DIMENSION(M,N) :: ZA,ZB,ZIN,ZOUT
      CALL DFFTW_INIT_THREADS()
      CALL DFFTW_PLAN_WITH_NTHREADS(PLAN,2)
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=CMPLX(DA)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZA=ZOUT
      ZIN=CMPLX(DB)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZB=ZOUT
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZA*CONJG(ZB)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DC=DBLE(ZOUT)
      RETURN
      END SUBROUTINE DCORR2D

C ******************************************************************************
      SUBROUTINE SPECKLEINTERFEROMETRY(TFILE,TRNG,RFILE,RRNG,PX,PY,
     &  PREFIX,BUFFERSIZE)
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER, INTENT(IN) :: PX,PY,TRNG(2),RRNG(2)
      INTEGER, OPTIONAL, INTENT(IN) :: BUFFERSIZE
      INTEGER :: X,Y,NAXES(3)
      DOUBLE PRECISION, DIMENSION(PX,PY) :: DTPSD,DRPSD,DPSD,DACF
      DOUBLE COMPLEX :: ZIN(PX,PY),ZOUT(PX,PY)
      CHARACTER*(*), INTENT(IN) :: OBSFILE,REFFILE,PREFIX
      INTERFACE
      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      INTEGER, INTENT(OUT) :: NAXES(3)
      END SUBROUTINE IMAGESIZE
      SUBROUTINE GETPSD(FILENAME,FPIXELS,LPIXELS,DX,DY,DPSD,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3),DX,DY
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      DOUBLE PRECISION, INTENT(OUT) :: DPSD(DX,DY)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE GETPSD
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: ARRAY(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
      SUBROUTINE DFFTSHIFT(W,H,DX)
      INTEGER, INTENT(IN) :: W,H
      DOUBLE PRECISION, INTENT(INOUT) :: DX(W,H)
      END SUBROUTINE DFFTSHIFT
      END INTERFACE
C  Statements:
C  ===========
      CALL IMAGESIZE(TFILE,NAXES)
      IF(TRNG(2).LT.NAXES(3))THEN
        PRINT *,'request is out of range.'
        RETURN
      END IF
      IF(PRESENT(BUFFERSIZE))THEN
        CALL GETPSD(TFILE,(/1,1,TRNG(1)/),(/NAXES(1),NAXES(2),TRNG(2)/),
     &    PX,PY,DTPSD,BUFFERSIZE)
      ELSE
        CALL GETPSD(TFILE,(/1,1,TRNG(1)/),(/NAXES(1),NAXES(2),TRNG(2)/),
     &    PX,PY,DTPSD)
      END IF
      CALL WRITEIMAGE(TRIM(PREFIX)//'_tgabs.fits',
     &  (/1,1,1/),(/PX,PY,1/),DSQRT(DTPSD))
      PRINT *,'spectral amplitude of target: '//
     &  TRIM(PREFIX)//'_tgabs.fits'
      CALL IMAGESIZE(RFILE,NAXES)
      IF(RRNG(2).LT.NAXES(3))THEN
        PRINT *,'request is out of range.'
        RETURN
      END IF
      IF(PRESENT(BUFFERSIZE))THEN
        CALL GETPSD(RFILE,(/1,1,RRNG(1)/),(/NAXES(1),NAXES(2),RRNG(2)/),
     &    PX,PY,DRPSD,BUFFERSIZE)
      ELSE
        CALL GETPSD(RFILE,(/1,1,RRNG(1)/),(/NAXES(1),NAXES(2),RRNG(2)/),
     &    PX,PY,DRPSD)
      END IF
      CALL WRITEIMAGE(TRIM(PREFIX)//'_rfabs.fits',
     &  (/1,1,1/),(/PX,PY,1/),DSQRT(DRPSD))
      PRINT *,'spectral amplitude of reference: '//
     &  TRIM(PREFIX)//'_rfabs.fits'
      DO X=1,PX
        DO Y=1,PY
          IF (DABS(DRPSD(X,Y)) .GT. 0)THEN
            DPSD(X,Y)=DTPSD(X,Y)/DRPSD(X,Y)
          ELSE
            DPSD(X,Y)=DBLE(0)
          END IF
        END DO
      END DO
      CALL WRITEIMAGE(TRIM(PREFIX)//'_dmpsd.fits',
     &  (/1,1,1/),(/PX,PY,1/),DPSD)
      PRINT *,'demodulated power spectral density: '//
     &  TRIM(PREFIX)//'_dmpsd.fits'
      CALL DFFTW_PLAN_DFT_2D(PLAN,PX,PY,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=CMPLX(DPSD)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DACF=DREAL(ZOUT/DBLE(NAXES(1))/DBLE(NAXES(2)))
      CALL DFFTSHIFT(PX,PY,DACF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_acf.fits',
     &  (/1,1,1/),(/PX,PY,1/),DACF)
      PRINT *,'reconstructed auto-correlation function: '//
     &  TRIM(PREFIX)//'_acf.fits'
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE SPECKLEINTERFEROMETRY
C ******************************************************************************
      SUBROUTINE GETPSD(FILENAME,FPIXELS,LPIXELS,DX,DY,DPSD,BUFFERSIZE)
C  Purpose:
C  ========
C  Get the mean power spectral density of all frames in the given FITS file.
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN):: FPIXELS(3),LPIXELS(3),DX,DY
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      INTEGER :: W,H,INFO,K,NPIXELS,NFRAMES,LBUF,NBUF,L,L1,L2
      INTEGER*8 :: PLAN
      DOUBLE PRECISION, INTENT(OUT) :: DPSD(DX,DY)
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:)
      DOUBLE COMPLEX, DIMENSION(DX,DY) :: ZIN,ZOUT
      CHARACTER*(*), INTENT(IN) :: FILENAME
      INTERFACE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      END INTERFACE
C  Statements:
C  ===========
      W=LPIXELS(1)-FPIXELS(1)+1
      H=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=W*H
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      IF(PRESENT(BUFFERSIZE))THEN
        LBUF=INT(FLOOR(DBLE(BUFFERSIZE*1024*1024)/DBLE(NPIXELS*8)))
      ELSE
        LBUF=INT(FLOOR(DBLE(1024*1024)/DBLE(NPIXELS)))
      END IF
      NBUF=INT(CEILING(DBLE(NFRAMES)/DBLE(LBUF)))
      ALLOCATE(DBUF(W,H,LBUF))
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
        RETURN
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_IMPORT_SYSTEM_WISDOM(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_IMPORT_SYSTEM_WISDOM failed.'
      END IF
C     PRINT *,'Start planning.'
      CALL DFFTW_PLAN_DFT_2D(PLAN,D1,D2,ZIN,ZOUT,-1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
C     PRINT *,'Finished planning.'
      DPSD=DBLE(0)
      DO K=1,NBUF
        L1=(K-1)*LBUF+1
        L2=MIN(NFRAMES,K*LBUF)
        CALL READIMAGE(FILENAME,(/FPIXELS(1),FPIXELS(2),L1/),
     &    (/LPIXELS(1),LPIXELS(2),L2/),DBUF)
        DO L=1,L2+1-L1
          ZIN=CMPLX(DBLE(0))
          ZIN(1:W,1:H)=DBUF(:,:,L)*DBLE(NPIXELS)/SUM(DBUF(:,:,L))
          CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
          DPSD=DPSD+DREAL(ZOUT*CONJG(ZOUT))
        END DO
      END DO
      DPSD=DPSD/DBLE(NFRAMES)
      DEALLOCATE(DBUF)
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE GETPSD
C ******************************************************************************
      SUBROUTINE DECONVCLEAN(M,N,DG,DF,DH,DBETA,MNUMIT)
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: X,Y,XC,YC,K,L,MNUMIT
      DOUBLE PRECISION, INTENT(IN) :: DG(M,N),DH(M,N),DBETA
      DOUBLE PRECISION, INTENT(OUT) :: DF(M,N)
      DOUBLE PRECISION :: DRES(M,N),DPSF(M,N),DSNR,DR
      XC=INT(CEILING(0.5*DBLE(N+1)))
      YC=INT(CEILING(0.5*DBLE(M+1)))
      DR=DBLE(0.5)*DBLE(MIN(M,N))
      X=MAXLOC(MAXVAL(DH,1),2)
      Y=MAXLOC(MAXVAL(DH,2),1)
      DPSF=EOSHIFT(EOSHIFT(DH,Y-YC,DBLE(0),1),X-XC,DBLE(0),2)
      DPSF=DPSF/SUM(DPSF)
      DRES=DG
      L=1
      DO K=1,MNUMIT
        L=L+1
        X=MAXLOC(MAXVAL(DRES,1),2)
        Y=MAXLOC(MAXVAL(DRES,2),1)
        DF(Y,X)=DF(Y,X)+DRES(Y,X)*DBETA
        DRES=DRES-DRES(Y,X)*DBETA*EOSHIFT(EOSHIFT(DPSF,YC-Y,DBLE(0),1),
     &    XC-X,DBLE(0),2)
        IF (SUM(DRES) .LE. 0)THEN
          EXIT
        END IF
        IF (L .GT. 100)THEN
          L=1
          PRINT *,'Loop: ',K
          CALL GETSNR(M,N,DR,DSNR,DRES,'P4')
          WRITE (*,'(A,ES10.3)') ' Residual SNR:',DSNR
          WRITE (*,'(A,ES10.3,A,ES10.3,A,ES10.3)')
     &     ' Residual sum:',SUM(DRES),
     &      ', CLEAN sum:',SUM(DF),', Total sum:',SUM(DF)+SUM(DRES)
          WRITE (*,'(A,I3,A,I3,A,ES10.3,ES10.3)')
     &      ' Maxima: (',X,',',Y,'),',MAXVAL(DRES),DRES(Y,X)
        END IF
      END DO
      CALL WRITEIMAGE('CLEAN_RESIDUAL.FITS',(/1,1,1/),(/M,N,1/),DRES)
      RETURN
      END SUBROUTINE DECONVCLEAN
C ******************************************************************************
      SUBROUTINE APPENDIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: STATUS,UNIT,BLOCKSIZE,BITPIX,NAXIS,GROUP,RWMODE,EXISTS
      INTEGER :: NAXES(3)
      DOUBLE PRECISION, INTENT(IN) :: ARRAY(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      BLOCKSIZE=1
      BITPIX=-64
      NAXIS=3
      RWMODE=1
      CALL FTEXIST(FILENAME,EXISTS,STATUS)
      IF (EXISTS .EQ. 1)THEN
C       PRINT *,TRIM(FILENAME),' exists.'
        CALL FTIOPN(UNIT,FILENAME,RWMODE,STATUS)
        CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
        IF (NFOUND .NE. 3)THEN
          PRINT *,'APPENDIMAGE failed to read the NAXIS keywords.'
          RETURN
        ENDIF
        NAXES(3)=NAXES(3)+LPIXELS(3)-FPIXELS(3)+1
        CALL FTUKYJ(UNIT,'NAXIS3',NAXES(3),'length of data axis 3',
     &    STATUS)
      ELSE
C       PRINT *,TRIM(FILENAME),' does not exist.'
        CALL FTINIT(UNIT,FILENAME,BLOCKSIZE,STATUS)
        NAXES=LPIXELS-FPIXELS+1
        CALL FTPHPS(UNIT,BITPIX,NAXIS,NAXES,STATUS)
      END IF
      GROUP=1
      CALL FTPSSD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,ARRAY,STATUS)
      CALL FTCLOS(UNIT,STATUS)
      CALL FTFIOU(UNIT,STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE APPENDIMAGE

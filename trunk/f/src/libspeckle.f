      SUBROUTINE PRINTERROR(STATUS)
      INTEGER :: STATUS
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
      ENDSUBROUTINE PRINTERROR
C ******************************************************************************
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      INTEGER :: STATUS,UNIT,BLOCKSIZE
      CHARACTER*(*) :: FILENAME
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
      RETURN
      ENDSUBROUTINE DELETEFILE
C ******************************************************************************
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER :: STATUS,UNIT,BLOCKSIZE,BITPIX,NAXIS,GROUP
      INTEGER :: FPIXELS(*),LPIXELS(*),NAXES(3)
      DOUBLE PRECISION :: ARRAY(*)
      CHARACTER*(*) :: FILENAME
      STATUS=0
      CALL DELETEFILE(FILENAME,STATUS)
      CALL FTGIOU(UNIT,STATUS)
      BLOCKSIZE=1
      BITPIX=-64
      NAXIS=3
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      CALL FTINIT(UNIT,FILENAME,BLOCKSIZE,STATUS)
      CALL FTPHPS(UNIT,BITPIX,NAXIS,NAXES,STATUS)
      GROUP=1
      CALL FTPSSD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,ARRAY,STATUS)
      CALL FTCLOS(UNIT,STATUS)
      CALL FTFIOU(UNIT,STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      ENDSUBROUTINE WRITEIMAGE
C ******************************************************************************
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER :: STATUS,UNIT,READWRITE,BLOCKSIZE,NAXIS,NFOUND
      INTEGER :: GROUP,FPIXELS(*),LPIXELS(*),INCS(3),NAXES(3)
      DOUBLE PRECISION :: ARRAY(*),NULLVAL
      LOGICAL :: ANYF
      CHARACTER*(*) :: FILENAME
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      READWRITE=0
      CALL FTOPEN(UNIT,FILENAME,READWRITE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'READIMAGE FAILED TO READ THE NAXISN KEYWORDS.'
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
      ENDSUBROUTINE READIMAGE
C ******************************************************************************
      SUBROUTINE AVERAGE(FILENAME,FPIXELS,LPIXELS,ARRAY,BUFFER)
      INTEGER :: FPIXELS(*),LPIXELS(*),K,NPIXELS,F(3),L(3)
      CHARACTER*(*) :: FILENAME
      DOUBLE PRECISION :: ARRAY(*),BUFFER(*),DA,DSUM
      INTEGER :: NAXES(3)
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      F(1)=FPIXELS(1)
      F(2)=FPIXELS(2)
      L(1)=LPIXELS(1)
      L(2)=LPIXELS(2)
      DO K=1,NPIXELS
        ARRAY(K)=0
      ENDDO
      DA=1.0
      DO K=FPIXELS(3),LPIXELS(3)
        F(3)=K
        L(3)=K
        CALL READIMAGE(FILENAME,F,L,BUFFER)
        DSUM=SUM(BUFFER(1:NPIXELS))
        CALL DSCAL(NPIXELS,DBLE(NPIXELS)/DSUM,BUFFER,1)
        CALL DAXPY(NPIXELS,DA,BUFFER,1,ARRAY,1)
      ENDDO
      CALL DSCAL(NPIXELS,1/DBLE(NAXES(3)),ARRAY,1)
      RETURN
      ENDSUBROUTINE AVERAGE
C ******************************************************************************
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER*(*) :: PATH,BASENAME,EXTNAME
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
      ENDSUBROUTINE RESOLVEPATH
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
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB
      INTEGER :: M,N,I,J,K
      DOUBLE PRECISION :: X,Y,XC,YC,D
      DOUBLE PRECISION :: IMG(M,*),BG(M,*),B(*)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:)
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
      ENDSUBROUTINE BGFIT2P2
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
      DOUBLE PRECISION :: IMG(M,N),BG(M,N),X,Y,XC,YC,D,B(*)
      INTEGER :: M,N,I,J,K
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:),WORK(:)
      INTEGER :: NPARAMS,NSAMPLES,LWORK,INFO,LDA,LDB
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
      ENDSUBROUTINE BGFIT2P4
C ******************************************************************************
      SUBROUTINE GETSNR(M,N,D,SNR,IMG,FTMETHOD)
      INTEGER :: M,N,FTMETHOD,I,J,K,NPIXS,NSPLS,INFO
      DOUBLE PRECISION :: SNR,D,X,Y,XC,YC,PS,PN
      DOUBLE PRECISION :: IMG(M,N),CLN(M,N),WORK(M,N)
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:)
      DOUBLE COMPLEX :: FFT2IN(M,N)
      DOUBLE COMPLEX :: COMM2D(M*N+3*(M+N)+100)
      DOUBLE COMPLEX, ALLOCATABLE :: FFTIN(:),COMM1D(:)
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
      ALLOCATE(FFTIN(NSPLS))
      ALLOCATE(COMM1D(3*NSPLS+100))
      IF (FTMETHOD .EQ. 2) THEN
        CALL BGFIT2P2(M,N,D,IMG,WORK,BUFFER)
      ELSE IF (FTMETHOD .EQ. 4) THEN
        CALL BGFIT2P4(M,N,D,IMG,WORK,BUFFER)
      ELSE IF (FTMETHOD .EQ. 0) THEN
        CALL BGFIT2P0(M,N,D,IMG,WORK)
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
      FFT2IN(1:M,1:N)=CMPLX(CLN(1:M,1:N))
      CALL ZFFT2D(-1,M,N,FFT2IN,COMM2D,INFO)
      PS=DZNRM2(M*N,RESHAPE(FFT2IN,(/M*N,1/)),1)
      PS=PS*PS/DBLE(NPIXS)
      FFTIN(1:NSPLS)=CMPLX(BUFFER(1:NSPLS))
      CALL ZFFT1D(-1,NSPLS,FFTIN,COMM1D,INFO)
      PN=DZNRM2(NSPLS,FFTIN,1)
      PN=PN*PN/DBLE(NSPLS)
      SNR=PS/PN
      DEALLOCATE(FFTIN)
      DEALLOCATE(BUFFER)
      DEALLOCATE(COMM1D)
      RETURN
      END SUBROUTINE GETSNR
C
C
C ******************************************************************************
C
C
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
      INTEGER :: M,N,INFO
      DOUBLE PRECISION :: DG(M,N),DF(M,N),DH(M,N)
      DOUBLE PRECISION :: DSNR
      DOUBLE COMPLEX :: ZG(M,N),ZF(M,N),ZH(M,N)
      DOUBLE COMPLEX :: ZCOMM(M*N+3*(M+N)),ZDECONV(M,N)
      ZG=CMPLX(DG)
      ZH=CMPLX(DH)
      CALL ZFFTSHIFT(M,N,ZH)
      CALL ZFFT2D(-1,M,N,ZG,ZCOMM,INFO)
      CALL ZFFT2D(-1,M,N,ZH,ZCOMM,INFO)
      ZDECONV=CONJG(ZH)/(ZH*CONJG(ZH)+CMPLX(DBLE(1)/DSNR))
      ZF=ZG*ZDECONV
      CALL ZFFT2D(1,M,N,ZF,ZCOMM,INFO)
      DF=DBLE(ZF)
      RETURN
      END SUBROUTINE DECONVWNR
C
C
C ******************************************************************************
C
C
      SUBROUTINE ZFFTSHIFT(M,N,ZX)
C  Shift zero-frequency component to centre of spectrum.
      INTEGER :: M,N,EM,SM
      DOUBLE COMPLEX :: ZX(M,N),ZY(M,N)
      SM=INT(FLOOR(REAL(M)*0.5))
      EM=INT(FLOOR(REAL(N)*0.5))
      ZY(1:M-SM,1:N-EM)=ZX(SM+1:M,EM+1:N)
      ZY(1:M-SM,N-EM+1:N)=ZX(SM+1:M,1:EM)
      ZY(M-SM+1:M,N-EM+1:N)=ZX(1:SM,1:EM)
      ZY(M-SM+1:M,1:N-EM)=ZX(1:SM,EM+1:N)
      ZX=ZY
      RETURN
      END SUBROUTINE ZFFTSHIFT
C
C
C ******************************************************************************
C
C
      SUBROUTINE SSAREC(INFILE,FPIXELS,LPIXELS,IMG,BG,WORK)
      INTEGER :: FPIXELS(*),LPIXELS(*),NAXES(2)
      INTEGER :: K,NPIXELS,NFRAMES,X,Y,XC,YC
      DOUBLE PRECISION :: IMG(*),BG(*),DA,DSUM
      DOUBLE PRECISION :: WORK(LPIXELS(1)-FPIXELS(1)+1,*)
      CHARACTER*(*) :: INFILE
      NAXES=(/LPIXELS(1)-FPIXELS(1)+1,LPIXELS(2)-FPIXELS(2)+1/)
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      XC=INT(CEILING(0.5*DBLE(NAXES(2)+1)))
      YC=INT(CEILING(0.5*DBLE(NAXES(1)+1)))
      CALL DSCAL(NPIXELS,DBLE(0),IMG,1)
      DO K=FPIXELS(3),LPIXELS(3)
          CALL READIMAGE(INFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &    (/LPIXELS(1),LPIXELS(2),K/),WORK)
        DSUM=SUM(WORK(1:NAXES(1),1:NAXES(2)))
        CALL DSCAL(NPIXELS,DBLE(NPIXELS)/DSUM,WORK,1)
        CALL DAXPY(NPIXELS,DBLE(-1),BG,1,WORK,1)
        X=MAXLOC(MAXVAL(WORK(1:NAXES(1),1:NAXES(2)),1),2)
        Y=MAXLOC(MAXVAL(WORK(1:NAXES(1),1:NAXES(2)),2),1)
        DSUM=SUM(WORK(1:NAXES(1),1:NAXES(2)))
        CALL DSCAL(NPIXELS,DBLE(NPIXELS)/DSUM,WORK,1)
        WORK(1:NAXES(1),1:NAXES(2))=EOSHIFT(
     &    WORK(1:NAXES(1),1:NAXES(2)),X-XC,DBLE(0),2)
        WORK(1:NAXES(1),1:NAXES(2))=EOSHIFT(
     &    WORK(1:NAXES(1),1:NAXES(2)),Y-YC,DBLE(0),1)
        CALL DAXPY(NPIXELS,DBLE(1),WORK,1,IMG,1)
      END DO
      CALL DSCAL(NPIXELS,DBLE(1)/DBLE(NFRAMES),IMG,1)
      RETURN
      END SUBROUTINE SSAREC
C
C
C ******************************************************************************
C
C
      SUBROUTINE ISAREC(OBSFILE,FPIXELS,LPIXELS,M,N,DREF,DBG,DISA,DR,
     &  P,NUMIT,PREFIX)
C
C  Arguments:
C  ==========
C  OBSFILE - Observed fits file to process.
C  FPIXELS - First pixels in each frame of the observed fits file.
C  LPIXELS - Last pixels in each frame of the observed fits file.
C  M       - Number of rows of DREF, DBG, and DISA.
C  N       - Number of columns of DREF, DBG, and DISA.
C  DREF    - Image of the reference star.
C  DBG     - Background of the average of all frames of the observed fits file.
C  DISA    - Input/Output. As input it is the initial estimate.
C  DR      - Radius of the border of the signal in the image.
C  P       - Pad all matrices to P x P for fast fourier transform.
C  NUMIT   - Number of iterations.
C  PREFIX  - Prefix of output filename.
C
C
      INTEGER :: FPIXELS(*),LPIXELS(*),NUMIT,K,L,INFO,X,Y
      INTEGER :: M,N,NPIXS,NUMFRM,P
      DOUBLE PRECISION :: DSNR,DR,DRMS
      DOUBLE PRECISION, DIMENSION(M,N) :: DIMG,DREF,DBG,DISA,DCORE
      DOUBLE PRECISION, DIMENSION(M,N) :: WORK
      DOUBLE PRECISION, DIMENSION(P,P) :: DCORR
      DOUBLE COMPLEX,DIMENSION(P,P) :: ZIMG,ZISA,ZPDS

      DOUBLE COMPLEX :: COMM(P*P+6*P+200)
      CHARACTER*(*) :: PREFIX
      CHARACTER :: OBSFILE*80,ISAFILE*80,NUMSTR*80
      NPIXS=M*N
      NUMFRM=LPIXELS(3)-FPIXELS(3)+1
      DO L=1,NUMIT
        WRITE(*,'(A,I4)') 'Loop =',L
        WRITE(NUMSTR,*) L
        CALL GETSNR(M,N,DR,DSNR,DISA,0)
        CALL DECONVWNR(M,N,DISA,DCORE,DREF,DSNR)
        ISAFILE=TRIM(PREFIX)//'_ISA_CORE_'//
     &    TRIM(ADJUSTL(NUMSTR))//'.FITS'
        CALL WRITEIMAGE(ISAFILE,(/1,1,1/),(/M,N,1/),DCORE)
        WORK=0.0
        ZISA=CMPLX(0.0)
        ZISA(1:M,1:N)=CMPLX(DCORE)
        CALL ZFFT2D(-1,P,P,ZISA,COMM,INFO)
        ZISA=CONJG(ZISA)
        DO K=FPIXELS(3),LPIXELS(3)
          CALL READIMAGE(OBSFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &      (/LPIXELS(1),LPIXELS(2),K/),DIMG)
          DIMG=DIMG/SUM(DIMG)*DBLE(NPIXS)-DBG
          DIMG=DIMG/SUM(DIMG)*DBLE(NPIXS)
          ZIMG=CMPLX(0.0)
          ZIMG(1:M,1:N)=CMPLX(DIMG)
          CALL ZFFT2D(-1,P,P,ZIMG,COMM,INFO)
          ZPDS=ZIMG*ZISA
          CALL ZFFT2D(1,P,P,ZPDS,COMM,INFO)
          DCORR=DBLE(ZPDS)
C         CALL WRITEIMAGE('CORR.FITS',(/1,1,1/),(/P,P,1/),DCORR)
          X=MAXLOC(MAXVAL(DCORR,1),2)
          Y=MAXLOC(MAXVAL(DCORR,2),1)
          IF (P-X.GT.X-1)THEN
            X=X-1
          ELSE
            X=X-P
          END IF
          IF (P-Y.GT.Y-1)THEN
            Y=Y-1
          ELSE
            Y=Y-P
          END IF
          WORK=WORK+EOSHIFT(EOSHIFT(DIMG,Y-1,DBLE(0),1),
     &      X-1,DBLE(0),2)
        END DO
        WORK=WORK/DBLE(NUMFRM)
        ISAFILE=TRIM(PREFIX)//'_ISA_'//TRIM(ADJUSTL(NUMSTR))//'.FITS'
        CALL WRITEIMAGE(ISAFILE,(/1,1,1/),(/M,N,1/),WORK)
        DRMS=SQRT(SUM((WORK-DISA)*(WORK-DISA))/DBLE(NPIXS))
        WRITE(*,'(A,A,ES10.3,A,ES10.3)') TRIM(ISAFILE),', RMS =',DRMS,
     &    ' SNR =',DSNR
        DISA=WORK
      END DO
      RETURN
      END SUBROUTINE
C
C
C ******************************************************************************
C
C
      SUBROUTINE DCORR2D(M,N,DA,DB,DC)
C  correlation coefficients of 2-dimensional array.
C
C  Purpose:
C  ========
C  c = corr(a, b). calculates c given a and b using fast fourier transform.
C  fft(c) = fft(a) * conj(fft(b)) (periodic extension is implied).
C
      INTEGER :: M,N,INFO
      DOUBLE PRECISION, DIMENSION(M,N) :: DA,DB,DC
      DOUBLE COMPLEX, DIMENSION(M,N) :: ZA,ZB,ZC
      DOUBLE COMPLEX :: COMM(M*N+3*(M+N))
      ZA=CMPLX(DA)
      ZB=CMPLX(DB)
      CALL ZFFT2D(-1,M,N,ZA,COMM,INFO)
      CALL ZFFT2D(-1,M,N,ZB,COMM,INFO)
      ZC=ZA*CONJG(ZB)
      CALL ZFFT2D(1,M,N,ZC,COMM,INFO)
      DC=DBLE(ZC)
      RETURN
      END SUBROUTINE DCORR2D
C
C
C ******************************************************************************
C
C
      SUBROUTINE BGFIT2P0(M,N,D,IMG,BG)
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB
      INTEGER :: M,N,I,J,K
      DOUBLE PRECISION :: X,Y,XC,YC,D,B
      DOUBLE PRECISION :: IMG(M,N),BG(M,N)
      NSAMPLES=0
      B=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
            B=IMG(I,J)
          END IF
        END DO
      END DO
      BG=B/DBLE(NSAMPLES)
      RETURN
      ENDSUBROUTINE BGFIT2P0

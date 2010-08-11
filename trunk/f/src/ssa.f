      PROGRAM MAIN
C  Simple shift-and-add routine.
C  =============================
C
C  Usage:
C  ======
C  ssa file_obs first_row,first_col,first_frm last_row,last_col,last_frm
C    radius prefix
C
C  Argument:
C  =========
C  1.    file_obs   - observed fits file to process.
C  2.1   first_row  - first row in each frame of file_obs.
C  2.2   first_col  - first column in each frame of file_obs.
C  2.3   first_frm  - first frame of file_obs to process.
C  3.1   last_row   - last row in each frame of file_obs.
C  3.2   last_col   - last column in each frame of file_obs.
C  3.3   last_frm   - last frame of file_obs to process.
C  4.    radius     - radius of the border of signal on the image.
C  5.    prefix     - prefix of output filename.
C
C  Declarations:
C  =============
      INTEGER :: STATUS,UNIT,FPIXELS(3),LPIXELS(3),NAXES(3)
      INTEGER :: K,NPIXELS
      DOUBLE PRECISION :: DR,DSNR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DAVG(:),
     &  DCLN(:),DB(:)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DIMG(:,:),
     &  DBG(:,:)
      CHARACTER*(256) :: INFILE,PREFIX,AVGFILE,ARG,BASENAME,EXTNAME,
     &  BGFILE,CLNFILE,OUTFILE
C  Statements:
C  ===========
      STATUS=0
C    Resolve the command line options:
C    =================================
      CALL GETARG(1,INFILE)
      CALL GETARG(2,ARG)
      READ(ARG,*) FPIXELS(1),FPIXELS(2),FPIXELS(3)
      CALL GETARG(3,ARG)
      READ(ARG,*) LPIXELS(1),LPIXELS(2),LPIXELS(3)
      CALL GETARG(4,ARG)
      READ(ARG,*) DR
      CALL GETARG(5,PREFIX)
C    Determine the size of the problem:
C    ==================================
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      PRINT *,'Definition:'
      PRINT *,'  Input: ',TRIM(INFILE)
      WRITE(*,'(A,I6)') '   First frame: ',FPIXELS(3)
      WRITE(*,'(A,I6)') '   Last frame: ',LPIXELS(3)
      WRITE(*,'(A,I5,A,I5,A)') '   First pixel: (',FPIXELS(1),',',
     &  FPIXELS(2),')'
      WRITE(*,'(A,I5,A,I5,A)') '   Last pixel: (',LPIXELS(1),',',
     &  LPIXELS(2),')'
      WRITE(*,'(A,I5,A,I5)') '   Frame size: ',NAXES(1),' x ',NAXES(2)
      WRITE(*,'(A,I6)') '   Number of frames: ',NAXES(3)
C    Allocate the space:
C    ===================
      ALLOCATE(DAVG(NPIXELS))
      ALLOCATE(DIMG(NAXES(1),NAXES(2)))
      ALLOCATE(DBG(NAXES(1),NAXES(2)))
      ALLOCATE(DCLN(NPIXELS))
      ALLOCATE(DB(NPIXELS))
C    Calculate the average of all the frames:
C    ========================================
      PRINT *,'Average of the frames:'
      CALL AVERAGE(INFILE,FPIXELS,LPIXELS,DAVG,100)
      AVGFILE=TRIM(PREFIX)//'_AVG.FITS'
      CALL WRITEIMAGE(AVGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DAVG)
      PRINT *,'  Output: ',TRIM(AVGFILE)
C    Fitting the background of the image:
C    ====================================
      PRINT *,'Calculate the background:'
      DIMG=RESHAPE(DAVG,(/NAXES(1),NAXES(2)/))
C      Using constant:
C      ======================
      PRINT *,'  Constant:'
      CALL BGFIT2P0(NAXES(1),NAXES(2),DR,DAVG,DBG,DB(1))
      PRINT *,'    Fitting result:'
      WRITE(*,'(A,ES10.3)') '     a_0 = ',DB(1)
      BGFILE=TRIM(PREFIX)//'_BG_P0.FITS'
      CALL WRITEIMAGE(BGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'    Background: ',TRIM(BGFILE)
      DCLN=DAVG(1:NPIXELS)-
     &  RESHAPE(DBG(1:NAXES(1),1:NAXES(2)),(/NPIXELS/))
      CLNFILE=TRIM(PREFIX)//'_CLN_P0.FITS'
      DCLN=DCLN/DBLE(NPIXELS)*SUM(DCLN)
      CALL WRITEIMAGE(CLNFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DCLN)
      PRINT *,'    Clean image: ',TRIM(CLNFILE)
      CALL GETSNR(NAXES(1),NAXES(2),DR,DSNR,DIMG,0)
      WRITE(*,'(A,ES10.3)') '     SNR: ',DSNR
C      Using 2nd polynomials:
C      ======================
      PRINT *,'  2nd Polynomials:'
      CALL BGFIT2P2(NAXES(1),NAXES(2),DR,DIMG,DBG,DB)
      PRINT *,'    Fitting result:'
      DO K=1,6
        WRITE(*,'(A,ES10.3,A,I1)') '     a_i = ',DB(K),', i = ',K
      END DO
      BGFILE=TRIM(PREFIX)//'_BG_P2.FITS'
      CALL WRITEIMAGE(BGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'    Background: ',TRIM(BGFILE)
      DCLN(1:NPIXELS)=DAVG(1:NPIXELS)-
     &  RESHAPE(DBG(1:NAXES(1),1:NAXES(2)),(/NPIXELS/))
      CLNFILE=TRIM(PREFIX)//'_CLN_P2.FITS'
      DCLN=DCLN/DBLE(NPIXELS)*SUM(DCLN)
      CALL WRITEIMAGE(CLNFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DCLN)
      PRINT *,'    Clean image: ',TRIM(CLNFILE)
      CALL GETSNR(NAXES(1),NAXES(2),DR,DSNR,DIMG,2)
      WRITE(*,'(A,ES10.3)') '     SNR: ',DSNR
C      Using 4th polynomials:
C      ======================
      PRINT *,'  4th Polynomials:'
      CALL BGFIT2P4(NAXES(1),NAXES(2),DR,DIMG,DBG,DB)
      PRINT *,'    Fitting result:'
      DO K=1,15
        WRITE(*,'(A,ES10.3,A,I2)') '     a_i = ',DB(K),', i = ',K
      END DO
      BGFILE=TRIM(PREFIX)//'_BG_P4.FITS'
      CALL WRITEIMAGE(BGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'    Background: ',TRIM(BGFILE)
      DCLN=DAVG(1:NPIXELS)-
     &  RESHAPE(DBG(1:NAXES(1),1:NAXES(2)),(/NPIXELS/))
      CLNFILE=TRIM(PREFIX)//'_CLN_P4.FITS'
      DCLN=DCLN/DBLE(NPIXELS)*SUM(DCLN)
      CALL WRITEIMAGE(CLNFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DCLN)
      PRINT *,'    Clean image: ',TRIM(CLNFILE)
      CALL GETSNR(NAXES(1),NAXES(2),DR,DSNR,DIMG,4)
      WRITE(*,'(A,ES10.3)') '     SNR: ',DSNR
C    Simple shift-and-add:
C    =====================
      PRINT *,'Simple shift-and-add:'
      CALL SSAREC(INFILE,FPIXELS,LPIXELS,DIMG,DBG)
      OUTFILE=TRIM(PREFIX)//'_SSA.FITS'
      CALL WRITEIMAGE(OUTFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      PRINT *,'    Output: ',TRIM(OUTFILE)
C    Deallocate space:
C    =================
      DEALLOCATE(DIMG)
      DEALLOCATE(DCLN)
      DEALLOCATE(DBG)
      DEALLOCATE(DAVG)
      DEALLOCATE(DB)
      STOP
      END PROGRAM MAIN

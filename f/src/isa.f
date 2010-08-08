      PROGRAM MAIN
C  Iterative shift-and-add routine.
C  ================================
C
C  Usage:
C  ======
C  isa file_obs first_row,first_col,first_frm last_row,last_col,last_frm
C    file_ref file_bg file_init radius pad_size num_it prefix
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
C  4.    file_ref   - pre-processed reference star fits file.
C  5.    file_bg    - estimate background fits file.
C  6.    file_init  - initial estimate fits file.
C  7.    radius     - radius of the border of signal on the image.
C  8.    pad_size   - pad matrices to the size p x p for fast fourier transform.
C  9.    num_it     - number of iterations.
C 10.    prefix     - prefix of output filename.
C
C  Declarations:
C  =============
      INTEGER :: FPIXELS(3),LPIXELS(3),M,N,NUMIT,P
      DOUBLE PRECISION :: DR
      DOUBLE PRECISION, ALLOCATABLE :: DREF(:,:),DISA(:,:),DBG(:,:)
      CHARACTER*(256) :: OBSFILE,REFFILE,INITFILE,PREFIX,ARG,BGFILE
C
C  Resolve the command line arguments:
C  ===================================
      CALL GETARG(1,OBSFILE)
      CALL GETARG(2,ARG)
      READ(ARG,*) FPIXELS(1),FPIXELS(2),FPIXELS(3)
      CALL GETARG(3,ARG)
      READ(ARG,*) LPIXELS(1),LPIXELS(2),LPIXELS(3)
      CALL GETARG(4,REFFILE)
      CALL GETARG(5,BGFILE)
      CALL GETARG(6,INITFILE)
      CALL GETARG(7,ARG)
      READ(ARG,*) DR
      CALL GETARG(8,ARG)
      READ(ARG,*) P
      CALL GETARG(9,ARG)
      READ(ARG,*) NUMIT
      CALL GETARG(10,PREFIX)
C
C  Main routine starts here:
C  =========================
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      ALLOCATE(DREF(M,N))
      ALLOCATE(DBG(M,N))
      ALLOCATE(DISA(M,N))
      CALL READIMAGE(REFFILE,(/1,1,1/),(/M,N,1/),DREF)
      CALL READIMAGE(BGFILE,(/1,1,1/),(/M,N,1/),DBG)
      CALL READIMAGE(INITFILE,(/1,1,1/),(/M,N,1/),DISA)
      CALL ISAREC(OBSFILE,FPIXELS,LPIXELS,M,N,DREF,DBG,DISA,DR,
     &  P,NUMIT,PREFIX)
      DEALLOCATE(DREF)
      DEALLOCATE(DBG)
      DEALLOCATE(DISA)
      STOP
      END PROGRAM MAIN

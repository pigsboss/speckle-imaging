      PROGRAM MAIN
C Speckle interferometry.
C =======================
C
C Usage:
C ======
C  si file_obs first_row_obs,first_col_obs,first_frm_obs
C    last_row_obs,last_col_obs,last_frm_obs radius_obs
C    file_ref first_row_ref,first_col_ref,first_frm_ref
C    last_row_ref,last_col_ref,last_frm_ref radius_ref
C    pad_size fit_method prefix
C
C Arguments:
C ==========
C  1.  file_obs      - observed file (target).
C  2.1 first_row_obs - first row number.
C  2.2 first_col_obs - first column number.
C  2.3 first_frm_obs - first frame number.
C  3.1 last_row_obs  - last row number.
C  3.2 last_col_obs  - last column number.
C  3.3 last_frm_obs  - last frame number.
C  4.  radius_obs    - radius of the border of signal in observed frames.
C  5.  file_ref      - reference file (reference star).
C  6.1 first_row_ref - first row number.
C  6.2 first_col_ref - first column number.
C  6.3 first_frm_ref - first frame number.
C  7.1 last_row_ref  - last row number.
C  7.2 last_col_ref  - last column number.
C  7.3 last_frm_ref  - last frame number.
C  8.  radius_ref    - radius of the border of signal in reference frames.
C  9.  pad_size      - pad all frames to this size before executing FFT.
C 10.  fit_method    - P0, P2, or P4.
C 11.  prefix        - prefix of output filename.
C
      INTEGER :: FPOBS(3),LPOBS(3),FPREF(3),LPREF(3),P,FTMETHOD
      CHARACTER(LEN=256) :: OBSFILE,REFFILE,PREFIX,ARG
      DOUBLE PRECISION :: DROBS,DRREF
      DOUBLE PRECISION, ALLOCATABLE :: DACF(:,:)
      CALL GETARG(1,OBSFILE)
      CALL GETARG(2,ARG)
      READ(ARG,*) FPOBS(1),FPOBS(2),FPOBS(3)
      CALL GETARG(3,ARG)
      READ(ARG,*) LPOBS(1),LPOBS(2),LPOBS(3)
      CALL GETARG(4,ARG)
      READ(ARG,*) DROBS
      CALL GETARG(5,REFFILE)
      CALL GETARG(6,ARG)
      READ(ARG,*) FPREF(1),FPREF(2),FPREF(3)
      CALL GETARG(7,ARG)
      READ(ARG,*) LPREF(1),LPREF(2),LPREF(3)
      CALL GETARG(8,ARG)
      READ(ARG,*) DRREF
      CALL GETARG(9,ARG)
      READ(ARG,*) P
      CALL GETARG(10,ARG)
      IF (TRIM(ARG) .EQ. 'P0')THEN
        FTMETHOD=0
      ELSE IF (TRIM(ARG) .EQ. 'P2')THEN
        FTMETHOD=2
      ELSE IF (TRIM(ARG) .EQ. 'P4')THEN
        FTMETHOD=4
      ELSE
        PRINT *,'Unknown fitting method ',TRIM(ARG)
        RETURN
      END IF
      CALL GETARG(11,PREFIX)
      ALLOCATE(DACF(P,P))
      CALL SIREC(OBSFILE,FPOBS,LPOBS,DROBS,
     &  REFFILE,FPREF,LPREF,DRREF,P,FTMETHOD,DACF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_ACF.FITS',
     &  (/1,1,1/),(/P,P,1/),DACF)
      DEALLOCATE(DACF)
      STOP
      END PROGRAM MAIN
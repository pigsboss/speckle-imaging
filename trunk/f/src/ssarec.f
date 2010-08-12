      SUBROUTINE SSAREC(INFILE,FPIXELS,LPIXELS,IMG,BG)
      INTEGER, INTENT(IN) :: FPIXELS(*),LPIXELS(*)
      INTEGER :: K,NPIXELS,NFRAMES,X,Y,XC,YC,NAXES(2)
      DOUBLE PRECISION, INTENT(OUT) :: IMG(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(2)-FPIXELS(2)+1)
      DOUBLE PRECISION, INTENT(IN) :: BG(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(2)-FPIXELS(2)+1)
      DOUBLE PRECISION :: WORK(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(2)-FPIXELS(2)+1),DB
      CHARACTER*(*), INTENT(IN) :: INFILE
      NAXES=(/LPIXELS(1)-FPIXELS(1)+1,LPIXELS(2)-FPIXELS(2)+1/)
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      XC=INT(CEILING(0.5*DBLE(NAXES(2)+1)))
      YC=INT(CEILING(0.5*DBLE(NAXES(1)+1)))
      IMG=DBLE(0)
      DO K=FPIXELS(3),LPIXELS(3)
          CALL READIMAGE(INFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &    (/LPIXELS(1),LPIXELS(2),K/),WORK)
        WORK=WORK/SUM(WORK)*DBLE(NPIXELS)-BG
        X=MAXLOC(MAXVAL(WORK,1),2)
        Y=MAXLOC(MAXVAL(WORK,2),1)
        WORK=WORK/SUM(WORK)*DBLE(NPIXELS)
        WORK=EOSHIFT(WORK,X-XC,DBLE(0),2)
        WORK=EOSHIFT(WORK,Y-YC,DBLE(0),1)
        IMG=IMG+WORK
      END DO
      IMG=IMG/DBLE(NFRAMES)
      RETURN
      END SUBROUTINE SSAREC
      PROGRAM MAIN
      INTEGER :: STATUS,UNIT,FPIXELS(3),LPIXELS(3),NAXES(3),NPIXELS
      CHARACTER :: FILENAME*80,ARG*255
      DOUBLE PRECISION, ALLOCATABLE :: ARRAY(:)
      STATUS=0
      CALL GETARG(1,FILENAME)
      WRITE(*,*) 'INPUT: ',FILENAME
      CALL GETARG(2,ARG)
      READ(ARG,*) FPIXELS(1)
      CALL GETARG(3,ARG)
      READ(ARG,*) FPIXELS(2)
      CALL GETARG(4,ARG)
      READ(ARG,*) FPIXELS(3)
      CALL GETARG(5,ARG)
      READ(ARG,*) LPIXELS(1)
      CALL GETARG(6,ARG)
      READ(ARG,*) LPIXELS(2)
      CALL GETARG(7,ARG)
      READ(ARG,*) LPIXELS(3)
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      PRINT *,'FRAME SIZE: ',NAXES(1),' BY ',NAXES(2)
      PRINT *,'NUMBER OF FRAMES: ',NAXES(3)
      ALLOCATE(ARRAY(NPIXELS))

C  Start the Simple shift-and-add routine.
      PRINT *,'START AVERAGING...'
      CALL AVERAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      PRINT *,'DONE.'
      CALL GETARG(8,FILENAME)
      WRITE(*,*) 'OUTPUT: ',FILENAME
      FPIXELS(1)=1
      FPIXELS(2)=1
      FPIXELS(3)=1
      LPIXELS(1)=NAXES(1)
      LPIXELS(2)=NAXES(2)
      LPIXELS(3)=1
      CALL WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      DEALLOCATE(ARRAY)
      STOP
      END
C **********************************************************************
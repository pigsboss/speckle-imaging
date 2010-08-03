      PROGRAM MAIN
      DOUBLE PRECISION A
      CHARACTER ARG*255
      CALL GETARG(1,ARG)
      READ(ARG,*) A
      WRITE(*,*) DSQRT(A)
      STOP
      END
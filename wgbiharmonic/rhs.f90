      REAL(8) FUNCTION rhs(x,y,testproblem)

      IMPLICIT NONE

      real(8) :: x,y
      real(8) :: pi = 3.14159265358979323846
      integer :: testproblem

      GO TO(100,200,300,400,500,600), testproblem

100   CONTINUE
      rhs = 8
      GOTO 900
      
200   CONTINUE
      rhs=4.d0*pi**4*DSIN(pi*x)*DSIN(pi*y)
      GOTO 900
      
300   CONTINUE
      rhs=4.d0*pi**4*DSIN(pi*x)*DCOS(pi*y)    
      GOTO 900

400   CONTINUE
      rhs=8.d0*pi**4*(3.d0-5.d0*DCOS(pi*y)**2-5.d0*DCOS(pi*x)**2+ &
          8.d0*DCOS(pi*x)**2*DCOS(pi*y)**2)
      GOTO 900

500   CONTINUE
      rhs=8.d0-48.d0*x-48.d0*y-48.d0*y**3+24.d0*y**4-48.d0*x**3+ &
          24.d0*x**4+288.d0*x**2*y**2+288.d0*x*y-288.d0*x*y**2-  &
          288.d0*x**2*y+72.d0*x**2+72.d0*y**2
      GOTO 900

600   CONTINUE
      rhs = 0.0

900   CONTINUE
      RETURN
      END FUNCTION rhs

      REAL(8) FUNCTION Exact_u(x,y,testproblem)

      IMPLICIT NONE

      real(8) :: x,y,u
      integer :: testproblem
      integer :: ix,iy
      real(8) :: pi=3.14159265358979323
      

      GO TO (10,20,30,40,50,60,70), testproblem
10    CONTINUE
      Exact_u=x*(1.d0-x)*y*(1-y)
      GOTO 90

20    CONTINUE
      Exact_u=DSIN(pi*x)*DSIN(pi*y)
      GOTO 90

30    CONTINUE
      Exact_u=DSIN(pi*x)*DCOS(pi*y)
      GOTO 90

40    CONTINUE
      u=DSIN(pi*x)*DSIN(pi*y)
      Exact_u=u*u
      GOTO 90

50    CONTINUE
      u=x*(1.d0-x)*y*(1.d0-y)
      Exact_u=u*u
      GOTO 90

60    CONTINUE
      Exact_u=((x**2+y**2)**(3.d0/4.d0))*(DSIN((3.d0/2.d0)* &
              DATAN(y/x))-3.d0*DSIN((1.d0/2.d0)*DATAN(y/x)))
      IF(x==0.0.AND.y==0.0) Exact_u = 0.d0
      GOTO 90
70    CONTINUE
      Exact_u=x*(1.d0-x)*y*(1.d0-y)

90    CONTINUE


      RETURN
      END FUNCTION Exact_u

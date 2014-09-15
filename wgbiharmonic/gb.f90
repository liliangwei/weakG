!------------------------------------------------
      REAL(8) FUNCTION Gb(x,y,testproblem,normal1,normal2)
      
      IMPLICIT NONE
      
      real(8) :: x,y,gux,guy,normal1,normal2
      integer :: testproblem
      real(8) :: pi = 3.1415926535897932384

      GO TO(10,20,30,40,50,60), testproblem
10    CONTINUE
      gux=(1.d0-2.d0*x)*y*(1.d0-y)
      guy=(1.d0-2.d0*y)*x*(1.d0-x)
      GOTO 90

20    CONTINUE
      gux=pi*DCOS(pi*x)*DSIN(pi*y)
      guy=pi*DSIN(pi*x)*DCOS(pi*y)
      GOTO 90

30    CONTINUE
      gux= pi*DCOS(pi*x)*DCOS(pi*y)
      guy=-pi*DSIN(pi*x)*DSIN(pi*y)
      GOTO 90

40    CONTINUE
      gux=2.d0*DSIN(pi*x)*DSIN(pi*y)**2*DCOS(pi*x)*pi
      guy=2.d0*DSIN(pi*x)**2*DSIN(pi*y)*DCOS(pi*y)*pi
      GOTO 90

50    CONTINUE
      gux=2.d0*x*(1.d0-x)**2*y**2*(1.d0-y)**2- &
          2.d0*x**2*(1.d0-x)*y**2*(1.d0-y)**2

      guy=2.d0*x**2*(1.d0-x)**2*y*(1.d0-y)**2 - &
          2.d0*x**2*(1.d0-x)**2*y**2*(1.d0-y)
      GOTO 90
60    CONTINUE
      gux=(3.d0/2.d0)*(x*DSIN((3.d0/2.d0)*DATAN(y/x))- &
          3.d0*x*DSIN((1.d0/2.d0)*DATAN(y/x))- &
          DCOS((3.d0/2.d0)*DATAN(y/x))*y+DCOS((1.d0/2.d0)* &
          DATAN(y/x))*y)/((x**2+y**2)**(1.d0/4.d0))
     
      guy=(3.d0/2.d0)*(x*DCOS((3.d0/2.d0)*DATAN(y/x))- &
          x*DCOS((1.d0/2.d0)*DATAN(y/x))+ &
          DSIN((3.d0/2.d0)*DATAN(y/x))*y-3.d0*DSIN((1.d0/2.d0)* &
          DATAN(y/x))*y)/((x**2+y**2)**(1.d0/4.d0))

90    CONTINUE

      Gb = gux*normal1+guy*normal2
      RETURN
      END FUNCTION Gb
!------------------------------------------------

      SUBROUTINE PLOT(u,ur)

      USE SHAREDATA
      IMPLICIT NONE
      character (len=90) :: filename
      real(8) :: u(n_node),ur(n_node)

      write(filename, '(A6,I1,A4)') 'Result',iter,'.dat'
      OPEN(100,FILE=filename,STATUS='REPLACE')

      write(100,*) 'TITLE = "FEWG 2D TRIANGULER MESH"'
      write(100,*) 'VARIABLES = "x", "y","u_cal","u_ture"'
      write(100,*) 'ZONE T="REFINED MESHES"   N=', n_node,  &
                   "E=",n_elem,"   F=FEPOINT", "   ET=TRIANGLE" 

      DO 110 i=1,n_node
      write(100,*) node(i,1),node(i,2),u(i),ur(i)
110   CONTINUE

      DO 120 i=1,n_elem
      write(100,*) (elem(i,k),k=1,3)
120   CONTINUE

      CLOSE(100)

      RETURN
      END

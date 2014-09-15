!------------------------------------------------
!   Weak Gulerkin method on biharmonic equation
!   From Dr.Lin Mu Matlab source code
!   Liangwei Li and Chunlei Liang
!         Jul/12/2014
!   
!   ' This is the first adaptation of WG method'
!------------------------------------------------

      PROGRAM WGRT02D

      USE SHAREDATA

!-----Divide [0,1]x[0,1] domain to 2xnxn triagles
!---- node = zeros((n+1)^2,2)
!---- elem = zeros(2*n^2,3)
!-------Control value determine the testcase-----
      testproblem = 7
!------------------------------------------------
      iter = 0

      DO 1000 nn=2,5

      iter = iter + 1

      n = 2**nn

      DO 1020 i=1,n+1
      DO 1021 j=1,n+1

!------Denote the coordinates

      node((i-1)*(n+1)+j,2)=1.0-(i-1)*1.0/n
      node((i-1)*(n+1)+j,1)=0.0+(j-1)*1.0/n

1021  CONTINUE
1020  CONTINUE
      n_node = n*(n+1)+n+1   ! size of node

      DO 1030 i=1,n
      DO 1031 j=1,n  

!--------triangulation connectivity
      elem((i-1)*2*n+2*(j-1)+1,1) = (i-1)*(n+1)+j
      elem((i-1)*2*n+2*(j-1)+1,2) = (i)*(n+1)+j
      elem((i-1)*2*n+2*(j-1)+1,3) = (i-1)*(n+1)+j+1
      elem((i-1)*2*n+2*(j-1)+2,1) = (i)*(n+1)+j+1
      elem((i-1)*2*n+2*(j-1)+2,2) = (i-1)*(n+1)+j+1
      elem((i-1)*2*n+2*(j-1)+2,3) = (i)*(n+1)+j

1031  CONTINUE
1030  CONTINUE

      n_elem = (n-1)*2*n+2*(n-1)+2

!-------Derive Mesh Information------------------
      CALL MESH

!-------Solve the equation-----------------------

!   Assemble stiffness matrix
      CALL STIFF

! --- Error Calculation
      Newerror(iter) = uerror
      L2error(iter)=uL2

      hh(iter) = EMAX

      write(*,*) hh(iter),Newerror(iter),L2error(iter)
      write(*,*) 'iteration = ',iter+1
      

1000  CONTINUE


      STOP
      END PROGRAM WGRT02D

!------------------------------------------------


      SUBROUTINE STIFF

      USE SHAREDATA

      IMPLICIT NONE

      integer              :: NT,NE,NO,dof,n_free  !Element, Edge, Node
      real(8), allocatable :: AA(:,:),area(:)
      real(8), allocatable :: aaa(:,:),bbb(:,:),ccc(:,:)
      real(8), allocatable :: x1(:),x2(:),x3(:)
      real(8), allocatable :: y1(:),y2(:),y3(:)
      real(8), allocatable :: l1(:),l2(:),l3(:)
      real(8), allocatable :: l12(:),l23(:),l13(:),l123(:)
      integer, allocatable :: n1(:),n2(:),n3(:),T(:)

      integer, allocatable :: Iv(:,:),Iu(:,:)
      real(8), allocatable :: val(:,:),b(:)
      real(8), allocatable :: ExactU(:),NDmid(:,:)
      real(8), allocatable :: tempAA(:,:),tempx(:)
      real(8)              :: tempxmid,tempymid
!
      real(8), allocatable :: x(:)
      integer, allocatable :: free(:)
      real(8)              :: Exact_u         ! Functions
      logical, allocatable :: isBdary(:)
      real(8), allocatable :: freeAA(:,:),invfreeAA(:,:)
      real(8), allocatable :: freex(:),freeb(:)

      real(8), allocatable :: xx(:),yy(:),w(:)
      real(8), allocatable :: quad_xy(:,:),quad_w(:)
      real(8), allocatable :: detAA(:,:)
      real(8)              :: rhs,det
      

      NT = n_elem                     ! Number of elements
      NE = n_edge                     ! Number of edges
      NO = n_node                     ! Number of nodes

      dof = NT+NE                     ! Degree of freedom

      allocate (AA(dof,dof))
      allocate (n1(NT),n2(NT),n3(NT))
      allocate (x1(NT),x2(NT),x3(NT))
      allocate (y1(NT),y2(NT),y3(NT))
      allocate (aaa(NT,3),bbb(NT,3),ccc(NT,3))
      allocate (l1(NT),l2(NT),l3(NT))
      allocate (l12(NT),l23(NT),l13(NT),l123(NT))
      allocate (area(NT))
!      allocate (E1(NT),E2(NT),E3(NT))

      allocate (Iv(NT,9),Iu(NT,9),val(NT,9))
      allocate (b(dof),ExactU(dof))

      DO 3111 i=1,NT
      n1(i) = elem(i,1)        
      n2(i) = elem(i,2)
      n3(i) = elem(i,3)
3111  CONTINUE

! ---  Coordinate of nodes
      DO 3112 i=1,NT
      x1(i) = node(n1(i),1)
      x2(i) = node(n2(i),1)
      x3(i) = node(n3(i),1)

      y1(i) = node(n1(i),2)
      y2(i) = node(n2(i),2)
      y3(i) = node(n3(i),2)
3112  CONTINUE
      
! --   Length of edges
      DO 3000 i=1,NT
      aaa(i,1) = y2(i)-y3(i)
      aaa(i,2) = y3(i)-y1(i)
      aaa(i,3) = y1(i)-y2(i)

      bbb(i,1) = x3(i)-x2(i)
      bbb(i,2) = x1(i)-x3(i)
      bbb(i,3) = x2(i)-x1(i)

      ccc(i,1) = x2(i)*y3(i) - x3(i)*y2(i)
      ccc(i,2) = x3(i)*y1(i) - x1(i)*y3(i)
      ccc(i,3) = x1(i)*y2(i) - x2(i)*y1(i)

      l1(i) = aaa(i,1)**2+bbb(i,1)**2
      l2(i) = aaa(i,2)**2+bbb(i,2)**2
      l3(i) = aaa(i,3)**2+bbb(i,3)**2

      area(i) = (-x2(i)*y1(i)+x3(i)*y1(i)+x1(i)*y2(i)  &
                -x3(i)*y2(i)-x1(i)*y3(i)+x2(i)*y3(i))/2.0d0

3000  CONTINUE

      l12 = l1 + l2
      l13 = l1 + l3
      l23 = l2 + l3
      l123 = l1 + l2 + l3

      allocate (T(NT))

      DO i=1,NT
      T(i) = i
      ENDDO
!
!  ---  Construct the matrix-----
!
      AA = 0.d0
!  ---Add the Matrix M00
      DO 4000 i=1,NT
      AA(T(i),T(i)) = AA(T(i),T(i)) + 144.d0*area(i)/l123(i)
4000  CONTINUE

!  ---Add the Matrix M0b and Mb0
      DO 4100 i=1,NT
      DO 4100 j=1,3
      Iu(i,j) = NT+ele2edge(T(i),j)
      Iv(i,j) = T(i)
      val(i,j) = -48.d0*area(i)/l123(i)
4100  CONTINUE

      DO 4200 i=1,NT
      DO 4200 j=1,3
      AA(Iu(i,j),Iv(i,j)) = AA(Iu(i,j),Iv(i,j))+val(i,j)
      AA(Iv(i,j),Iu(i,j)) = AA(Iv(i,j),Iu(i,j))+val(i,j)
4200  CONTINUE

! ---Add the Matrix Mbb

      Iu = 0.d0; Iv = 0.d0; val = 0.d0;
      DO 4300 j=1,9
      DO 4300 i=1,NT
      IF(j<=3) Iu(i,j) = NT+ele2edge(T(i),1)
      IF(j>=4.AND.j<=6) Iu(i,j) = NT+ele2edge(T(i),2)
      IF(j>=7.AND.j<=9) Iu(i,j) = NT+ele2edge(T(i),3)

      IF(MOD(j,3)==1) Iv(i,j) = NT+ele2edge(T(i),1)
      IF(MOD(j,3)==2) Iv(i,j) = NT+ele2edge(T(i),2)
      IF(MOD(j,3)==0) Iv(i,j) = NT+ele2edge(T(i),3)
4300  CONTINUE

      DO 4400 i=1,NT
      val(i,1)=16.d0*area(i)/l123(i)+2.d0*l1(i)/(2.d0*area(i))
      val(i,2)=16.d0*area(i)/l123(i)+(l3(i)-l12(i))/(2.d0*area(i))
      val(i,3)=16.d0*area(i)/l123(i)+(l2(i)-l13(i))/(2.d0*area(i))
      val(i,4)=16.d0*area(i)/l123(i)+(l3(i)-l12(i))/(2.d0*area(i))
      val(i,5)=16.d0*area(i)/l123(i)+2.d0*l2(i)/(2.d0*area(i))
      val(i,6)=16.d0*area(i)/l123(i)+(l1(i)-l23(i))/(2.d0*area(i))
      val(i,7)=16.d0*area(i)/l123(i)+(l2(i)-l13(i))/(2.d0*area(i))
      val(i,8)=16.d0*area(i)/l123(i)+(l1(i)-l23(i))/(2.d0*area(i))
      val(i,9)=16.d0*area(i)/l123(i)+2.d0*l3(i)/(2.d0*area(i))
4400  CONTINUE

      DO 4500 i=1,NT
      DO 4500 j=1,9
      AA(Iu(i,j),Iv(i,j)) = AA(Iu(i,j),Iv(i,j))+val(i,j)
4500  CONTINUE


!------------------------------------------------
! --- Assemble the RHS
! --- Volume force by Gaussian Integral
      quad_num = 13  ! denote the type
!------------------------------------------------

      allocate(quad_xy(quad_num,2),quad_w(quad_num))

      CALL QUAD_RULE(quad_num,quad_xy,quad_w)
! 
! --- Handle the Dirichlet Boundary condition----
!
      allocate (x(dof),isBdary(dof),free(dof))
      allocate (xx(quad_num),yy(quad_num),w(quad_num))
      allocate (detAA(3,3))
      allocate (NDmid(n_ex,2))

      detAA = 0.d0

      DO 4550 j=1,NT
      b(j) = 0.d0
      DO 4560 i=1,quad_num
      xx(i)=node(elem(j,1),1)*quad_xy(i,1)+node(elem(j,2),1)*  &
            quad_xy(i,2)+node(elem(j,3),1)*(1.d0-quad_xy(i,1)- &
            quad_xy(i,2))
      
      yy(i)=node(elem(j,1),2)*quad_xy(i,1)+node(elem(j,2),2)*  &
            quad_xy(i,2)+node(elem(j,3),2)*(1-quad_xy(i,1)-    &
            quad_xy(i,2))
      
      DO 4570 k=1,3
      detAA(1,k) = 1.d0
      detAA(2,k) = node(elem(j,k),1)
      detAA(3,k) = node(elem(j,k),2)
4570  CONTINUE 

      det=detAA(1,1)*(detAA(2,2)*detAA(3,3)-detAA(3,2)*detAA(2,3)) &
         +detAA(2,1)*(detAA(1,3)*detAA(3,2)-detAA(1,2)*detAA(3,3)) &
         +detAA(3,1)*(detAA(1,2)*detAA(2,3)-detAA(1,3)*detAA(2,2))

      w(i)=0.5d0*det*quad_w(i)
      b(j)=b(j)+w(i)*rhs(xx(i),yy(i),testproblem)
4560  CONTINUE
4550  CONTINUE
      
!---  Handle boundary conditions
      x = 0.d0;k=0;
      isBdary = .FALSE.   ! Initialize as all false
      DO 4575 i=1,n_ex
      isBdary(NT+exterioredge(i,4)) = .TRUE.
4575  CONTINUE

      DO 4580 i=1,n_ex
      DO 4580 j=1,2
      NDmid(i,j) = (node(exterioredge(i,1),j)+   &
                    node(exterioredge(i,2),j))/2.d0
4580  CONTINUE

      DO 4590 i=1,dof
      IF(isBdary(i)==.TRUE.) THEN
      k=k+1
      x(i) = Exact_u(NDmid(k,1),NDmid(k,2),testproblem)
      ENDIF
4590  CONTINUE

      k=0
      DO 4600 i=1,dof
      IF(isBdary(i)==.FALSE.) THEN
      k=k+1
      free(k) = i
      ENDIF
4600  CONTINUE
      n_free = k
      
      allocate(tempAA(dof,dof),tempx(dof))
      tempAA=0.d0;tempx=0.d0;k=0;
      DO 4610 i=1,dof
      IF(isBdary(i)==.TRUE.) THEN
      k=k+1
      DO 4620 j=1,dof
      tempx(k) = x(i)
      tempAA(j,k) = AA(j,i)
      b(j) = b(j) - tempAA(j,k)*tempx(k)
4620  CONTINUE
      ENDIF
4610  CONTINUE


! --- x(free = AA(free,free)\b(free)
      allocate(freeAA(n_free,n_free),freex(n_free),freeb(n_free))

      DO 4630 i=1,n_free
      freeb(i) = b(free(i))
      DO 4630 j=1,n_free
      freeAA(i,j) = AA(free(i),free(j))
4630  CONTINUE
       
      CALL GAUSS(n_free,freeAA,freeb,n_free,freex)

      DO 4650 i=1,n_free
      x(free(i)) = freex(i)
4650  CONTINUE

      ExactU = 0.d0
      DO 4660 i=1,NT
      tempxmid=1.d0/3.d0*(node(elem(i,1),1)+node(elem(i,2),1)+   &
                       node(elem(i,3),1))
      tempymid=1.d0/3.d0*(node(elem(i,1),2)+node(elem(i,2),2)+   &
                       node(elem(i,3),2))
      ExactU(i) = Exact_u(tempxmid,tempymid,testproblem)
4660  CONTINUE

      DO 4680 i=NT+1,NT+NE
      k = i-NT
      tempxmid=0.5d0*(node(edge(k,1),1)+node(edge(k,2),1))
      tempymid=0.5d0*(node(edge(k,1),2)+node(edge(k,2),2))
      ExactU(i) = Exact_u(tempxmid,tempymid,testproblem)
4680  CONTINUE


! --- Error report
      uerror = 0.d0; uL2=0.d0;
      DO 5000 i=1,dof
      DO 5000 j=1,dof
      uerror=uerror+(x(i)-ExactU(i))*AA(i,j)*(x(j)-ExactU(j))
5000  CONTINUE
      uerror = DSQRT(uerror)

      DO 5100 i=1,NT
      uL2=uL2+area(i)*(x(i)-ExactU(i))**2
5100  CONTINUE
      uL2 = DSQRT(uL2)


! --- Plot out
      CALL PLOT(x)
      
      RETURN
      END

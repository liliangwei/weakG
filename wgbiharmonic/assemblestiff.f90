
      SUBROUTINE STIFF

      USE SHAREDATA

      IMPLICIT NONE

      integer              :: NT,NE,NO  !Element, Edge, Node
      real(8), allocatable :: AA(:,:),area(:)
      real(8), allocatable :: aaa(:,:),bbb(:,:),ccc(:,:)
      real(8), allocatable :: E1(:),E2(:),E3(:)
      real(8), allocatable :: G1_1(:,:),G1_2(:,:),G2_1(:,:)
      real(8), allocatable :: G2_2(:,:),G3_1(:,:),G3_2(:,:)
      real(8), allocatable :: x1(:),x2(:),x3(:)
      real(8), allocatable :: y1(:),y2(:),y3(:)
      integer, allocatable          :: n1(:),n2(:),n3(:)

      integer, allocatable :: T(:)
      integer, allocatable :: it1(:),it2(:),it3(:)
      integer, allocatable :: itt1(:),itt2(:),itt3(:)
      integer, allocatable :: flag1(:),flag2(:),flag3(:)
      integer              :: temp1,temp2

      integer, allocatable :: Iv(:,:),Iu(:,:)
      integer, allocatable :: IvIu(:,:),IuIv(:,:)
      real(8), allocatable :: val(:,:),valval(:,:)
      real(8), allocatable :: b(:),Mid(:,:)  ! b matrix in matlab
      real(8), allocatable :: ExactU(:)
!
! --- boundary condition
      real(8), allocatable :: x(:),mid1(:,:),mid2(:,:)
      real(8), allocatable :: NDnode(:,:),NDmid(:,:)
      real(8), allocatable :: normal_el(:,:),normal(:,:)
      integer, allocatable :: free(:)
      integer              :: n_free
      integer, allocatable :: exterioredge4(:)
      real(8)              :: Exact_u,Gb    ! Functions
      logical, allocatable :: isBdary(:),isBdary1(:),isBdary2(:)
      logical, allocatable :: isBdary3(:),isBdary4(:)
      real(8), allocatable :: tempAA(:,:),tempx(:)
      real(8), allocatable :: freeAA(:,:),invfreeAA(:,:)
      real(8), allocatable :: freex(:),freeb(:)

      real(8), allocatable :: G1(:,:),G2(:,:)
      real(8), allocatable :: xx(:),yy(:),w(:)
      real(8), allocatable :: nabla1(:),nabla2(:),nabla3(:)
      real(8), allocatable :: uh(:),u(:),guhx(:),guhy(:)
      real(8)              :: tmp,tmpp
      real(8), allocatable :: gux(:),guy(:)
      real(8), allocatable :: p1(:),p2(:),p3(:),p4(:)
      real(8), allocatable :: p5(:),p6(:)
      real(8)              :: temp_M(6,6),temp_b(6),invtemp_M(6,6)
      real(8)              :: uu(6),uuu(13)
      real(8), allocatable :: length(:)
      real(8), allocatable :: quad_xy(:,:),quad_w(:);
      integer, allocatable :: BNode(:)

      real(8)              :: rhs
      

      NT = n_elem                     ! Number of elements
      NE = n_edge                     ! Number of edges
      NO = n_node                     ! Number of nodes

      allocate (AA(NO+3*NE,NO+3*NE))
      allocate (n1(NT),n2(NT),n3(NT))
      allocate (x1(NT),x2(NT),x3(NT))
      allocate (y1(NT),y2(NT),y3(NT))
      allocate (aaa(NT,3),bbb(NT,3),ccc(NT,3))
      allocate (area(NT))
      allocate (E1(NT),E2(NT),E3(NT))
      allocate (G1_1(NT,2),G1_2(NT,2),G2_1(NT,2))
      allocate (G2_2(NT,2),G3_1(NT,2),G3_2(NT,2))

      allocate (Iv(NT,36),Iu(NT,36),val(NT,36))
      allocate (IvIu(NT,72),IuIv(NT,72),valval(NT,72))
      allocate (b(NO+3*NE),Mid(NT,2))
      allocate (ExactU(NO+3*NE))

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

      area(i) = (-x2(i)*y1(i)+x3(i)*y1(i)+x1(i)*y2(i)  &
                -x3(i)*y2(i)-x1(i)*y3(i)+x2(i)*y3(i))/2.0d0

      E1(i) = DSQRT(aaa(i,1)**2+bbb(i,1)**2)   ! length of edge
      E2(i) = DSQRT(aaa(i,2)**2+bbb(i,2)**2)   ! length of edge
      E3(i) = DSQRT(aaa(i,3)**2+bbb(i,3)**2)   ! length of edge
3000  CONTINUE

!   -- Compute the Gaussian 2-points on each edge
      DO 3100 i=1,NT
      DO 3100 j=1,2

      G1_1(i,j)=node(elem(i,2),j)+(3.d0-DSQRT(3.d0))/6.d0* &
                (node(elem(i,3),j)-node(elem(i,2),j))
      G1_2(i,j)=node(elem(i,2),j)+(3.d0+DSQRT(3.d0))/6.d0* &
                (node(elem(i,3),j)-node(elem(i,2),j))

      G2_1(i,j)=node(elem(i,3),j)+(3.d0-DSQRT(3.d0))/6.d0* &
                (node(elem(i,1),j)-node(elem(i,3),j))
      G2_2(i,j)=node(elem(i,3),j)+(3.d0+DSQRT(3.d0))/6.d0* &
                (node(elem(i,1),j)-node(elem(i,3),j))


      G3_1(i,j)=node(elem(i,1),j)+(3.d0-DSQRT(3.d0))/6.d0* &
                (node(elem(i,2),j)-node(elem(i,1),j))
      G3_2(i,j)=node(elem(i,1),j)+(3.d0+DSQRT(3.d0))/6.d0* &
                (node(elem(i,2),j)-node(elem(i,1),j))
3100  CONTINUE

      allocate (T(NT),it1(NT),it2(NT),it3(NT))
      allocate (itt1(NT),itt2(NT),itt3(NT))
      allocate (flag1(NT),flag2(NT),flag3(NT))

      DO i=1,NT
      T(i) = i
      ENDDO
 
      it1 = 0;it2 = 0;it3 = 0;
      itt1 = 0;itt2 = 0;itt3 = 0;
      flag1 = 0;flag2 = 0;flag3 = 0;

! --- Get the order number of the Gaussian points on each edge
! --- here, it~ means 2-Points
! --- indexedge

      DO 3200 i=1,NT

      it1(i) = 2*ele2edge(i,1)-1
      it2(i) = 2*ele2edge(i,2)-1
      it3(i) = 2*ele2edge(i,3)-1

      itt1(i) = it1(i)+1
      itt2(i) = it2(i)+1
      itt3(i) = it3(i)+1

      flag1(i) = 1
      flag2(i) = 1
      flag3(i) = 1

      IF(elem(i,2).NE.edge(ele2edge(i,1),1)) THEN
      temp1 = it1(i)
      temp2 = itt1(i)
      it1(i) = temp2
      itt1(i) = temp1

      flag1(i) = -1
      ENDIF

      IF(elem(i,3).NE.edge(ele2edge(i,2),1)) THEN
      temp1 = it2(i)
      temp2 = itt2(i)
      it2(i) = temp2
      itt2(i) = temp1

      flag2(i) = -1
      ENDIF

      IF(elem(i,1).NE.edge(ele2edge(i,3),1)) THEN
      temp1 = it3(i)
      temp2 = itt3(i)
      it3(i) = temp2
      itt3(i) = temp1

      flag3(i) = -1
      ENDIF
3200  CONTINUE
    
!
!  ---  Construct the matrix-----
!

! ---   Compute the (\delta_w U_h, \delta_w V)
      DO 3300 j=1,36
      DO 3310 i=1,NT

      IF(j<=6)  Iv(i,j) = NO+NE+it1(i) 
      IF(j>=7.AND.j<=12) Iv(i,j) = NO+NE+itt1(i)
      IF(j>=13.AND.j<=18) Iv(i,j) = NO+NE+it2(i)
      IF(j>=19.AND.j<=24) Iv(i,j) = NO+NE+itt2(i)
      IF(j>=25.AND.j<=30) Iv(i,j) = NO+NE+it3(i)
      IF(j>=31.AND.j<=36) Iv(i,j) = NO+NE+itt3(i)

3310  CONTINUE
3300  CONTINUE

      DO 3400 j=1,36
      DO 3410 i=1,NT

      IF(MOD(j,6)==1) Iu(i,j) = NO+NE+it1(i)
      IF(MOD(j,6)==2) Iu(i,j) = NO+NE+itt1(i)
      IF(MOD(j,6)==3) Iu(i,j) = NO+NE+it2(i)
      IF(MOD(j,6)==4) Iu(i,j) = NO+NE+itt2(i)
      IF(MOD(j,6)==5) Iu(i,j) = NO+NE+it3(i)
      IF(MOD(j,6)==0) Iu(i,j) = NO+NE+itt3(i)

3410  CONTINUE
3400  CONTINUE

      DO 3500 i=1,NT

      val(i,1) = E1(i)**2/(4.d0*area(i))
      val(i,2) = E1(i)**2/(4.d0*area(i))
      val(i,3) = flag1(i)*flag2(i)*E1(i)*E2(i)/(4.d0*area(i))
      val(i,4) = val(i,3)
      val(i,5) = flag1(i)*flag3(i)*E1(i)*E3(i)/(4.d0*area(i))
      val(i,6) = val(i,5)

      val(i,7) = val(i,7-6)
      val(i,8) = val(i,8-6)
      val(i,9) = val(i,9-6)
      val(i,10) = val(i,10-6)
      val(i,11) = val(i,11-6)
      val(i,12) = val(i,12-6)

      val(i,13) = flag1(i)*flag2(i)*E1(i)*E2(i)/(4.d0*area(i))
      val(i,14) = val(i,13)
      val(i,15) = E2(i)**2/(4.d0*area(i))
      val(i,16) = val(i,15)
      val(i,17) = flag2(i)*flag3(i)*E2(i)*E3(i)/(4.d0*area(i))
      val(i,18) = val(i,17)
      
      val(i,19) = val(i,19-6)
      val(i,20) = val(i,20-6)
      val(i,21) = val(i,21-6)
      val(i,22) = val(i,22-6)
      val(i,23) = val(i,23-6)
      val(i,24) = val(i,24-6)
      
      val(i,25) = flag1(i)*flag3(i)*E1(i)*E3(i)/(4.d0*area(i))
      val(i,26) = val(i,25)
      val(i,27) = flag2(i)*flag3(i)*E2(i)*E3(i)/(4.d0*area(i))
      val(i,28) = val(i,27)
      val(i,29) = E3(i)**2/(4.d0*area(i))
      val(i,30) = val(i,29)

      val(i,31) = val(i,31-6)
      val(i,32) = val(i,32-6)
      val(i,33) = val(i,33-6)
      val(i,34) = val(i,34-6)
      val(i,35) = val(i,35-6)
      val(i,36) = val(i,36-6)
3500  CONTINUE

! --- Deonte matrix AA

      AA = 0.0d0

      DO 3600 i=1,NT
      DO 3600 j=1,36
      AA(Iv(i,j),Iu(i,j)) = AA(Iv(i,j),Iu(i,j))+val(i,j)
3600  CONTINUE

! ---  clear val Iu Iv
      Iu = 0;Iv = 0;val = 0.0d0;
!
! --- Add penalized term h^{-1}<(\nabla u_0-u_nb n_e)\cdot n,
! --- \nabla v_0-v_nb n_e)\cdot n> to AA
!
!   h^{-1}<\nabla u_0\cdot n,\nabla v_0\cdot n>_{partial K}
!   h^{-1}<\nabla u_0\cdot n,\nabla v_0\cdot n>
      
      DO 3700 j=1,36
      DO 3710 i=1,NT

      IF(j<=6)  Iu(i,j) = elem(i,1) 
      IF(j>=7.AND.j<=12) Iu(i,j) = elem(i,2)
      IF(j>=13.AND.j<=18) Iu(i,j) = elem(i,3)
      IF(j>=19.AND.j<=24) Iu(i,j) = NO+ele2edge(i,1)
      IF(j>=25.AND.j<=30) Iu(i,j) = NO+ele2edge(i,2)
      IF(j>=31.AND.j<=36) Iu(i,j) = NO+ele2edge(i,3)

3710  CONTINUE
3700  CONTINUE

      DO 3720 j=1,36
      DO 3730 i=1,NT

      IF(MOD(j,6)==1) Iv(i,j) = elem(i,1)
      IF(MOD(j,6)==2) Iv(i,j) = elem(i,2)
      IF(MOD(j,6)==3) Iv(i,j) = elem(i,3)
      IF(MOD(j,6)==4) Iv(i,j) = NO+ele2edge(i,1)
      IF(MOD(j,6)==5) Iv(i,j) = NO+ele2edge(i,2)
      IF(MOD(j,6)==0) Iv(i,j) = NO+ele2edge(i,3)

3730  CONTINUE
3720  CONTINUE


      DO 3760 i=1,NT
! (1,1)
      val(i,1)=(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+7.d0/3.d0*   &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+7.d0/3.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (1,2)
      val(i,2)=-(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-   &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-1.d0/3.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (1,3)
      val(i,3)=-(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))- 1.d0/3.d0*  &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (1,4)
      val(i,4)=-2.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))- 2.d0*  &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+2.d0/3.d0*   &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+2.d0/3.d0*  &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (1,5)
      val(i,5)=-2.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+2.d0/3.d0*  &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+10.d0/3.d0*   &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+10.d0/3.d0*  &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (1,6)
      val(i,6)=-2.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+10.d0/3.d0*  &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+2.d0/3.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))+10.d0/3.d0*  &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (2,1)
      val(i,7)=-(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-   &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-1.d0/3.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (2,2)
      val(i,8)=7.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+7.d0/3.d0*   &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (2,3)
      val(i,9)=-1.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+7.d0/3.d0*   &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (2,4) 
      val(i,10)=2.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+10.d0/3.d0*  &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-2.0d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+10.d0/3.d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (2,5)
      val(i,11)=2.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-2.0d0*  &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-2.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+2.d0/3.d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (2,6) 
      val(i,12)=10.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-2.0d0*  &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+10.d0/3.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))+2.d0/3.d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (3,1)
      val(i,13)=-(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-1.d0/3.d0*   &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (3,2)
      val(i,14)=-1.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-1.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-   &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (3,3)
      val(i,15)=7.d0/3.d0*(aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+7.d0/3.d0*   &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+   &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (3,4)  
      val(i,16)=10.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+2.d0/3.d0*  &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+10.d0/3.d0*   &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-2.0d0*  &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (3,5)
      val(i,17)=10.d0/3.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+10.d0/3.d0*  &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+2.d0/3.d0*   &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E3(i)**2))-2.0d0*  &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (3,6)
      val(i,18)=2.d0/3.d0*(aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+2.d0/3.d0*  &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-2.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))-2.0d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (4,1)     
      val(i,19)=-2.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-2.d0*  &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+2.d0/3.d0*   &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+2.0d0/3.d0*  &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (4,2)
      val(i,20)=2.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+10.d0/3.d0*  &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-2.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+10.0d0/3.d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (4,3)
      val(i,21)=10.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+2.d0/3.d0*  &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+10.d0/3.d0*   &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-2.0d0*  &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (4,4)
      val(i,22)=16.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-16.d0/3.d0*  &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+16.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+16.d0/3.d0*  &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (4,5)
      val(i,23)=8.d0/3.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1) &
               -aaa(i,1)*aaa(i,1)-bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+8.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1) &
               -aaa(i,2)*aaa(i,2)-bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+8.d0/3.d0*   &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (4,6)
      val(i,24)=8.d0/3.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1) &
               -aaa(i,1)*aaa(i,1)-bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+8.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+8.d0/3.d0*   &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3) &
               -aaa(i,3)*aaa(i,3)-bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (5,1)
      val(i,25)=-2.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+2.d0/3.d0*  &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+10.d0/3.d0*   &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+10.d0/3.d0*  &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (5,2)
      val(i,26)=2.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-2.d0*  &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-2.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+2.d0/3.d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (5,3)
      val(i,27)=10.d0/3.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+10.d0/3.d0*  &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+2.d0/3.d0*   &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))* &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-2.d0*  &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))

! (5,4)
      val(i,28)=8.d0/3.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1) &
               -aaa(i,1)*aaa(i,1)-bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+8.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1) &
               -aaa(i,2)*aaa(i,2)-bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+8.d0/3.d0*   &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))

! (5,5)
      val(i,29)=16.d0/3.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+16.d0/3.d0*  &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-16.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+16.d0/3.d0*  &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (5,6)
      val(i,30)=8.d0/3.d0*(aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3) &
               -aaa(i,3)*aaa(i,3)-bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))+8.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3) &
               -aaa(i,2)*aaa(i,2)-bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+8.d0/3.d0*   &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))
! (6,1)
      val(i,31)=-2.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+10.d0/3.d0*  &
               (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+2.d0/3.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))+10.d0/3.d0*  &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (6,2)
      val(i,32)=10.d0/3.d0*(aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))-2.d0*  &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+10.d0/3.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))+2.d0/3.d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (6,3)
      val(i,33)=2.d0/3.d0*(aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+2.d0/3.d0*  &
               (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))-2.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))-2.d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (6,4)
      val(i,34)=8.d0/3.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1) &
               -aaa(i,1)*aaa(i,1)-bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+8.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+8.d0/3.d0*   &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3) &
               -aaa(i,3)*aaa(i,3)-bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))
! (6,5)
      val(i,35)=8.d0/3.d0*(aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3) &
               -aaa(i,3)*aaa(i,3)-bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))+8.d0/3.d0*   &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3) &
               -aaa(i,2)*aaa(i,2)-bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+8.d0/3.d0*   &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))
! (6,6)
      val(i,36)=16.d0/3.d0*(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))* &
               (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
               /((4.d0*area(i)**2)*(E1(i)**2))+16.d0/3.d0*  &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))* &
               (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
               /((4.d0*area(i)**2)*(E2(i)**2))+16.d0/3.d0*   &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))* &
               (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))-16.d0/3.d0*  &
               (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))* &
               (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
               /((4.d0*area(i)**2)*(E3(i)**2))

3760  CONTINUE

      DO 3800 i=1,NT
      DO 3800 j=1,36
      AA(Iv(i,j),Iu(i,j)) = AA(Iv(i,j),Iu(i,j))+val(i,j)
3800  CONTINUE

!
! ---  h^{-1}<(u_nb n_e)\cdot n,(v_nb n_e)\cdot n>
! Since u_nb is P1 on each edge
 
! --- clear Iv and val
      Iv = 0;val = 0.d0;

      DO 3900 i=1,NT
      Iv(i,1) = NO+NE+it1(i)
      Iv(i,2) = NO+NE+itt1(i)
      Iv(i,3) = NO+NE+it2(i)
      Iv(i,4) = NO+NE+itt2(i)
      Iv(i,5) = NO+NE+it3(i)
      Iv(i,6) = NO+NE+itt3(i)
      
3900  CONTINUE

      DO 3910 i=1,NT
      DO 3910 j=1,6
      val(i,j) = 0.5d0
3910  CONTINUE

      DO 3920 i=1,NT
      DO 3920 j=1,6
      AA(Iv(i,j),Iv(i,j)) = AA(Iv(i,j),Iv(i,j))+val(i,j)
3920  CONTINUE

! 
!  -h^{-1}<\nabla u_0\cdot n,(v_nb n_e)\cdot n>_{\partial K}
! & h^{-1}<\nabla v_0\cdot n,\nabla u\cdot n>_{\partial K}
!  \nabla u\cdot n={G1,G2)*flag on each Gaussian Points
!

! --- clear Iu,Iv, val
      Iu = 0; Iv = 0; val = 0.d0;

      DO 3930 j=1,36
      DO 3930 i=1,NT

      IF(j<=6)  Iu(i,j) = elem(i,1) 
      IF(j>=7.AND.j<=12) Iu(i,j) = elem(i,2)
      IF(j>=13.AND.j<=18) Iu(i,j) = elem(i,3)
      IF(j>=19.AND.j<=24) Iu(i,j) = NO+ele2edge(i,1)
      IF(j>=25.AND.j<=30) Iu(i,j) = NO+ele2edge(i,2)
      IF(j>=31.AND.j<=36) Iu(i,j) = NO+ele2edge(i,3)

      IF(MOD(j,6)==1) Iv(i,j) = NO+NE+it1(i)
      IF(MOD(j,6)==2) Iv(i,j) = NO+NE+itt1(i)
      IF(MOD(j,6)==3) Iv(i,j) = NO+NE+it2(i)
      IF(MOD(j,6)==4) Iv(i,j) = NO+NE+itt2(i)
      IF(MOD(j,6)==5) Iv(i,j) = NO+NE+it3(i)
      IF(MOD(j,6)==0) Iv(i,j) = NO+NE+itt3(i)

3930  CONTINUE

      DO 3950 i=1,NT

      val(i,1) = -(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
                 /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,2) = -(aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1)) &
                 /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,3) = (3.d0-2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))  &
                 /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,4) = (3.d0+2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))  &
                 /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,5) = (3.d0+2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))  &
                 /(2.d0*area(i)*E3(i))*flag3(i)
      val(i,6) = (3.d0-2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))  &
                 /(2.d0*area(i)*E3(i))*flag3(i)

      
      val(i,7) = (3.d0+2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))  &
                 /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,8) = (3.d0-2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,1)*aaa(i,2)+bbb(i,1)*bbb(i,2))  &
                 /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,9) = -(aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
                 /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,10) = -(aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2)) &
                 /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,11) = (3.d0-2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))  &
                 /(2.d0*area(i)*E3(i))*flag3(i)
      val(i,12) = (3.d0+2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))  &
                 /(2.d0*area(i)*E3(i))*flag3(i)

      
      val(i,13) = (3.d0-2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))  &
                 /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,14) = (3.d0+2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,1)*aaa(i,3)+bbb(i,1)*bbb(i,3))  &
                 /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,15) = (3.d0+2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))  &
                 /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,16) = (3.d0-2.d0*DSQRT(3.d0))/3.d0*  &
                 (aaa(i,2)*aaa(i,3)+bbb(i,2)*bbb(i,3))  &
                 /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,17) = -(aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
                 /(2.d0*area(i)*E3(i))*flag3(i)
      val(i,18) = -(aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3)) &
                 /(2.d0*area(i)*E3(i))*flag3(i)

      
      val(i,19) = 2.d0/3.d0*((3.d0-DSQRT(3.d0))*  &
                  (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))+ &
                  (3.d0+DSQRT(3.d0))* &
                  (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))) &
                  /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,20) = 2.d0/3.d0*((3.d0+DSQRT(3.d0))*  &
                  (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))+ &
                  (3.d0-DSQRT(3.d0))* &
                  (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))) &
                  /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,21) = 2.d0/3.d0*(3.d0+DSQRT(3.d0)) * &
                  (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))  &
                  /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,22) = 2.d0/3.d0*(3.d0-DSQRT(3.d0)) * &
                  (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))  &
                  /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,23) = 2.d0/3.d0*(3.d0-DSQRT(3.d0)) * &
                  (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))  &
                  /(2.d0*area(i)*E2(i))*flag3(i)
      val(i,24) = 2.d0/3.d0*(3.d0+DSQRT(3.d0)) * &
                  (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))  &
                  /(2.d0*area(i)*E2(i))*flag3(i)


      val(i,25) = 2.d0/3.d0*(3.d0-DSQRT(3.d0)) * &
                  (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))  &
                  /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,26) = 2.d0/3.d0*(3.d0+DSQRT(3.d0)) * &
                  (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))  &
                  /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,27) = 2.d0/3.d0*((3.d0+DSQRT(3.d0))*  &
                  (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))+ &
                  (3.d0-DSQRT(3.d0))* &
                  (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))) &
                  /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,28) = 2.d0/3.d0*((3.d0-DSQRT(3.d0))*  &
                  (aaa(i,2)*aaa(i,1)+bbb(i,2)*bbb(i,1))+ &
                  (3.d0+DSQRT(3.d0))* &
                  (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))) &
                  /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,29) = 2.d0/3.d0*(3.d0+DSQRT(3.d0)) * &
                  (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))  &
                  /(2.d0*area(i)*E3(i))*flag3(i)
      val(i,30) = 2.d0/3.d0*(3.d0-DSQRT(3.d0)) * &
                  (aaa(i,3)*aaa(i,3)+bbb(i,3)*bbb(i,3))  &
                  /(2.d0*area(i)*E3(i))*flag3(i)


      val(i,31) = 2.d0/3.d0*(3.d0+DSQRT(3.d0)) * &
                  (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))  &
                  /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,32) = 2.d0/3.d0*(3.d0-DSQRT(3.d0)) * &
                  (aaa(i,1)*aaa(i,1)+bbb(i,1)*bbb(i,1))  &
                  /(2.d0*area(i)*E1(i))*flag1(i)
      val(i,33) = 2.d0/3.d0*(3.d0-DSQRT(3.d0)) * &
                  (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))  &
                  /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,34) = 2.d0/3.d0*(3.d0+DSQRT(3.d0)) * &
                  (aaa(i,2)*aaa(i,2)+bbb(i,2)*bbb(i,2))  &
                  /(2.d0*area(i)*E2(i))*flag2(i)
      val(i,35) = 2.d0/3.d0*((3.d0-DSQRT(3.d0))*  &
                  (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))+ &
                  (3.d0+DSQRT(3.d0))* &
                  (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))) &
                  /(2.d0*area(i)*E3(i))*flag3(i)
      val(i,36) = 2.d0/3.d0*((3.d0+DSQRT(3.d0))*  &
                  (aaa(i,3)*aaa(i,1)+bbb(i,3)*bbb(i,1))+ &
                  (3.d0-DSQRT(3.d0))* &
                  (aaa(i,3)*aaa(i,2)+bbb(i,3)*bbb(i,2))) &
                  /(2.d0*area(i)*E3(i))*flag3(i)

3950  CONTINUE

      val = -0.5d0*val

      
      DO 3980 j=1,72
      DO 3980 i=1,NT

      IF(j<=36) THEN
      IuIv(i,j) = Iu(i,j)
      IvIu(i,j) = Iv(i,j)
      valval(i,j) = -val(i,j)
      ELSE
      IuIv(i,j) = Iv(i,j-36)
      IvIu(i,j) = Iu(i,j-36)
      valval(i,j) = -val(i,j-36)
      ENDIF
3980  CONTINUE

      DO 3990 i=1,NT
      DO 3990 j=1,72
      AA(IuIv(i,j),IvIu(i,j)) = AA(IuIv(i,j),IvIu(i,j))+valval(i,j)
3990  CONTINUE

      DO 4000 i=1,NT
      DO 4000 j=1,2
      Mid(i,j) = (node(elem(i,1),j)+node(elem(i,2),j) &
                 +node(elem(i,3),j))/3.d0
4000  CONTINUE

      Iv = 0; val = 0.d0; b = 0.d0;

      DO 4100 i=1,NT
      Iv(i,1) = elem(i,1)
      Iv(i,2) = elem(i,2)
      Iv(i,3) = elem(i,3)
      Iv(i,4) = NO+ele2edge(i,1)
      Iv(i,5) = NO+ele2edge(i,2)
      Iv(i,6) = NO+ele2edge(i,3)
4100  CONTINUE

      DO 4150 i=1,NT
      val(i,1)=rhs(node(elem(i,1),1),node(elem(i,1),2), &
               testproblem)*area(i)/20.d0-rhs(Mid(i,1),Mid(i,2),&
               testproblem)*area(i)/20.d0

      val(i,2)=rhs(node(elem(i,2),1),node(elem(i,2),2), &
               testproblem)*area(i)/20.d0-rhs(Mid(i,1),Mid(i,2),&
               testproblem)*area(i)/20.d0

      val(i,3)=rhs(node(elem(i,3),1),node(elem(i,3),2), &
               testproblem)*area(i)/20.d0-rhs(Mid(i,1),Mid(i,2),&
               testproblem)*area(i)/20.d0

      val(i,4)=rhs(0.5*node(elem(i,2),1)+0.5*node(elem(i,3),1), &
               0.5*node(elem(i,2),2)+0.5*node(elem(i,3),2),     &
               testproblem)*area(i)/15.d0*2.d0+                 &
               rhs(Mid(i,1),Mid(i,2),testproblem)*area(i)/5.d0

      val(i,5)=rhs(0.5*node(elem(i,1),1)+0.5*node(elem(i,3),1), &
               0.5*node(elem(i,1),2)+0.5*node(elem(i,3),2),     &
               testproblem)*area(i)/15.d0*2.d0+                 &
               rhs(Mid(i,1),Mid(i,2),testproblem)*area(i)/5.d0

      val(i,6)=rhs(0.5*node(elem(i,2),1)+0.5*node(elem(i,1),1), &
               0.5*node(elem(i,2),2)+0.5*node(elem(i,1),2),     &
               testproblem)*area(i)/15.d0*2.d0+                 &
               rhs(Mid(i,1),Mid(i,2),testproblem)*area(i)/5.d0
4150  CONTINUE

      DO 4180 i=1,NT
      DO 4180 j=1,6
      b(Iv(i,j)) = b(Iv(i,j))+val(i,j) 
4180  CONTINUE
!------------------------------------------------
      quad_num = 13  ! denote the type
!------------------------------------------------

      allocate(quad_xy(quad_num,2),quad_w(quad_num))

      CALL QUAD_RULE(quad_num,quad_xy,quad_w)

! 
! --- Handle the Dirichlet Boundary condition----
!
      allocate (x(NO+3*NE),isBdary(NO+3*NE))
      allocate (exterioredge4(n_ex))

      isBdary = .FALSE.   ! Initialize as all false

      exterioredge4 = exterioredge(:,1)

      allocate(BNode(n_ex))
      CALL Unique(exterioredge4,n_ex,BNode)

      DO 4200 i=1,size(BNode)
      isBdary(BNode(i)) = .TRUE.
4200  CONTINUE
      
      DO 4210 i=1,n_ex
      isBdary(NO+exterioredge(i,4)) = .TRUE.
      isBdary(NO+NE+2*exterioredge(i,4)-1) = .TRUE.
      isBdary(NO+NE+2*exterioredge(i,4)) = .TRUE.
4210  CONTINUE

      allocate(isBdary1(MAXVAL(BNode)),isBdary2(NO+NE))
      allocate(isBdary3(NO+3*NE-1),isBdary4(NO+3*NE))

      isBdary1 = .FALSE.
      isBdary2 = .FALSE.

      DO 4220 i=1,size(BNode)
      isBdary1(BNode(i)) = .TRUE.
4220  CONTINUE

      DO 4230 i=1,n_ex
      isBdary2(NO+exterioredge(i,4)) = .TRUE.
4230  CONTINUE

      allocate(NDnode(n_ex,2),NDmid(n_ex,2))

      NDnode = 0.d0;NDmid = 0.d0;

      DO 4250 i=1,n_ex
      DO 4250 j=1,2
      NDnode(i,j) = node(BNode(i),j)
      NDmid(i,j) = (node(exterioredge(i,1),j) + &
                    node(exterioredge(i,2),j))/2.d0
4250  CONTINUE

      x = 0.d0; k = 0;
      DO 4300 i=1,size(isBdary1)
      IF(isBdary1(i)==.TRUE.) THEN
      k=k+1
      x(i) = Exact_u(NDnode(k,1),NDnode(k,2),testproblem)
      ENDIF
4300  CONTINUE

      k = 0
      DO 4310 i=1,size(isBdary2)
      IF(isBdary2(i)==.TRUE.) THEN
      k=k+1
      x(i) = Exact_u(NDmid(k,1),NDmid(k,2),testproblem)
      ENDIF
4310  CONTINUE

      isBdary3 = .FALSE.
      isBdary4 = .FALSE.

      DO 4320 i=1,n_ex
      isBdary3(NO+NE+2*exterioredge(i,4)-1) = .TRUE.
      isBdary4(NO+NE+2*exterioredge(i,4)) = .TRUE.
4320  CONTINUE

      allocate(mid1(n_ex,2),mid2(n_ex,2))
      
      DO 4330 i=1,n_ex
      DO 4330 j=1,2
      mid1(i,j) = (3.d0-DSQRT(3.d0))/6.d0*(node(exterioredge(i,2),j) &
              -node(exterioredge(i,1),j))+node(exterioredge(i,1),j)
      mid2(i,j) = (3.d0+DSQRT(3.d0))/6.d0*(node(exterioredge(i,2),j) &
              -node(exterioredge(i,1),j))+node(exterioredge(i,1),j)
4330  CONTINUE

      allocate(normal_el(n_ex,2))

      DO 4340 i=1,n_ex
      normal_el(i,1)=(node(exterioredge(i,2),2)- &
                     (node(exterioredge(i,1),2)))/ &
                     DSQRT((node(exterioredge(i,2),1) -&
                     node(exterioredge(i,1),1))**2 + &
                     (node(exterioredge(i,2),2)- &
                     node(exterioredge(i,1),2))**2)

      normal_el(i,2)=(node(exterioredge(i,1),1)- &
                     (node(exterioredge(i,2),1)))/ &
                     DSQRT((node(exterioredge(i,2),1) -&
                     node(exterioredge(i,1),1))**2 + &
                     (node(exterioredge(i,2),2)- &
                     node(exterioredge(i,1),2))**2)
4340  CONTINUE

      k=0
      DO 4350 i=1,size(isBdary3)
      IF(isBdary3(i)==.TRUE.) THEN
      k=k+1
      x(i) = Gb(mid1(k,1),mid1(k,2),testproblem,normal_el(k,1), &
                      normal_el(k,2))
      ENDIF
4350  CONTINUE

      k=0
      DO 4360 i=1,size(isBdary4)
      IF(isBdary4(i)==.TRUE.) THEN
      k=k+1
      x(i) = Gb(mid2(k,1),mid2(k,2),testproblem,normal_el(k,1), &
                      normal_el(k,2))
      ENDIF
4360  CONTINUE

      allocate(free(NO+3*NE))

      k=0
      DO 4370 i=1,NO+3*NE
      IF(isBdary(i)==.FALSE.) THEN
      k=k+1
      free(k) = i
      ENDIF
4370  CONTINUE
      n_free = k       

      allocate(tempAA(NO+3*NE,NO+3*NE),tempx(NO+3*NE))
      tempAA = 0.d0; tempx = 0.d0;k=0;

! --- b = b - AA(:,isBdary)*(isBdary)
      DO 4380 i=1,NO+3*NE
      IF(isBdary(i)==.TRUE.) THEN
      k=k+1
      DO 4390 j=1,NO+3*NE
      tempx(k) = x(i)
      tempAA(j,k) = AA(j,i)
      b(j) = b(j) - tempAA(j,k)*tempx(k)
4390  CONTINUE
      ENDIF
4380  CONTINUE

! --- x(free)=AA(free,free)\b(free)
! --- inverse AA to invAA

      allocate(freeAA(n_free,n_free),invfreeAA(n_free,n_free))
      allocate(freex(n_free),freeb(n_free))

      DO 4395 i=1,n_free
      DO 4395 j=1,n_free
      freeAA(i,j) = AA(free(i),free(j))     ! Extract AA(free,free)
4395  CONTINUE
      
      DO 4396 i=1,n_free
      freeb(i) = b(free(i))
4396  CONTINUE

      CALL GAUSS(n_free,freeAA,freeb,n_free,freex)

      DO 4397 i=1,n_free
      x(free(i))=freex(i)
4397  CONTINUE

!      CALL inverse(freeAA,invfreeAA,n_free) ! inverse AA(free,free)
!      
!      DO 4400 j=1,n_free
!      DO 4400 i=1,n_free
!      x(free(i)) = x(free(i))+invfreeAA(i,j)*b(free(j))
!4400  CONTINUE

      
      ExactU=0.d0
      DO 4500 i=1,NO
      ExactU(i) = Exact_u(node(i,1),node(i,2),testproblem)
4500  CONTINUE

      k=0
      DO 4550 i=1,NE
      k=i+NO
      ExactU(k)=Exact_u((node(edge(i,1),1)+node(edge(i,2),1))*  &
                    0.5d0,(node(edge(i,1),2)+node(edge(i,2),2)) &
                    *0.5d0,testproblem)
4550  CONTINUE   

      allocate(G1(NE,2),G2(NE,2))

      DO 4560 i=1,NE
      DO 4560 j=1,2
      G1(i,j)=(3.d0-DSQRT(3.d0))/6.d0*(node(edge(i,2),j)- &
               node(edge(i,1),j))+node(edge(i,1),j)
      G2(i,j)=(3.d0+DSQRT(3.d0))/6.d0*(node(edge(i,2),j)- &
               node(edge(i,1),j))+node(edge(i,1),j)
4560  CONTINUE

      allocate(normal(NE,2))
      DO 4570 i=1,NE
      normal(i,1)=(node(edge(i,2),2)-node(edge(i,1),2))/ &
                  DSQRT((node(edge(i,2),2)-node(edge(i,1),2))**2+ &
                  (node(edge(i,2),1)-node(edge(i,1),1))**2)
      normal(i,2)=(node(edge(i,1),1)-node(edge(i,2),1))/ &
                  DSQRT((node(edge(i,2),2)-node(edge(i,1),2))**2+ &
                  (node(edge(i,2),1)-node(edge(i,1),1))**2)
4570  CONTINUE

      k=NO+NE
      DO 4580 i=1,NE
      ExactU(k+2*i-1)=Gb(G1(i,1),G1(i,2),testproblem,normal(i,1),&
                     normal(i,2))
      ExactU(k+2*i)=Gb(G2(i,1),G2(i,2),testproblem,normal(i,1),&
                     normal(i,2))
4580  CONTINUE
      
      uL2=0.d0; uH1=0.d0; uQL2=0.d0; u0Inf=0.d0;
      
      allocate(xx(quad_num),yy(quad_num),w(quad_num))
      allocate(nabla1(quad_num),nabla2(quad_num),nabla3(quad_num))
      allocate(uh(quad_num),u(quad_num))
      allocate(gux(quad_num),guy(quad_num))
      allocate(guhx(quad_num),guhy(quad_num))
      allocate(p1(quad_num),p2(quad_num),p3(quad_num))
      allocate(p4(quad_num),p5(quad_num),p6(quad_num))


      
      DO 5000 j=1,NT
      temp_M = 0.0; temp_b = 0.0; invtemp_M=0.0;uu=0.0;
      DO 5100 i=1,quad_num

      xx(i)=node(elem(j,1),1)*quad_xy(i,1)+node(elem(j,2),1)* &
            quad_xy(i,2)+node(elem(j,3),1)*(1-quad_xy(i,1)- &
            quad_xy(i,2))
      yy(i)=node(elem(j,1),2)*quad_xy(i,1)+node(elem(j,2),2)* &
            quad_xy(i,2)+node(elem(j,3),2)*(1-quad_xy(i,1)- &
            quad_xy(i,2)) 
      w(i)=0.5d0*area(j)*quad_w(i)

      nabla1(i)=(aaa(j,1)*xx(i)+bbb(j,1)*yy(i)+ccc(j,1))/ &
                 (2.d0*area(j))
      nabla2(i)=(aaa(j,2)*xx(i)+bbb(j,2)*yy(i)+ccc(j,2))/ &
                 (2.d0*area(j))
      nabla3(i)=(aaa(j,3)*xx(i)+bbb(j,3)*yy(i)+ccc(j,3))/ &
                 (2.d0*area(j))

      uh(i)=x(elem(j,1))*nabla1(i)*(2*nabla1(i)-1)+  &
            x(elem(j,2))*nabla2(i)*(2*nabla2(i)-1)+  &
            x(elem(j,3))*nabla3(i)*(2*nabla3(i)-1)+  &
            x(NO+ele2edge(j,1))*4.d0*nabla2(i)*nabla3(i)+ &
            x(NO+ele2edge(j,2))*4.d0*nabla1(i)*nabla3(i)+ &
            x(NO+ele2edge(j,3))*4.d0*nabla2(i)*nabla1(i)

      u(i)=Exact_u(xx(i),yy(i),testproblem)

      uL2=uL2+w(i)*(u(i)-uh(i))**2

      if(i==1) THEN 
      tmp=DABS(uh(i)-u(i))
      ELSE
      tmpp=DABS(uh(i)-u(i))
      IF(tmpp.GT.tmp) tmp=tmpp
      ENDIF

      u0Inf=max(tmp,u0Inf)

      CALL Gradient_u(xx(i),yy(i),testproblem,gux(i),guy(i))

      guhx(i)=x(elem(j,1))*(4.d0*nabla1(i)-1)*aaa(j,1)/  &
               (2.d0*area(j))+                          &
              x(elem(j,2))*(4.d0*nabla2(i)-1)*aaa(j,2)/  &
               (2.d0*area(j))+                          &
              x(elem(j,3))*(4.d0*nabla3(i)-1)*aaa(j,3)/  &
               (2.d0*area(j))+                          &
              x(NO+ele2edge(j,1))*(4.d0*nabla2(i)*aaa(j,3)+ &
              4.d0*nabla3(i)*aaa(j,2))/(2.d0*area(j))+      &
              x(NO+ele2edge(j,2))*(4.d0*nabla1(i)*aaa(j,3)+ &
              4.d0*nabla3(i)*aaa(j,1))/(2.d0*area(j))+       &
              x(NO+ele2edge(j,3))*(4.d0*nabla2(i)*aaa(j,1)+ &
              4.d0*nabla1(i)*aaa(j,2))/(2.d0*area(j))       
              
      guhy(i)=x(elem(j,1))*(4.d0*nabla1(i)-1)*bbb(j,1)/  &
               (2.d0*area(j))+                          &
              x(elem(j,2))*(4.d0*nabla2(i)-1)*bbb(j,2)/  &
               (2.d0*area(j))+                          &
              x(elem(j,3))*(4.d0*nabla3(i)-1)*bbb(j,3)/  &
               (2.d0*area(j))+                          &
              x(NO+ele2edge(j,1))*(4.d0*nabla2(i)*bbb(j,3)+ &
              4.d0*nabla3(i)*bbb(j,2))/(2.d0*area(j))+      &
              x(NO+ele2edge(j,2))*(4.d0*nabla1(i)*bbb(j,3)+ &
              4.d0*nabla3(i)*bbb(j,1))/(2.d0*area(j))+       &
              x(NO+ele2edge(j,3))*(4.d0*nabla2(i)*bbb(j,1)+ &
              4.d0*nabla1(i)*bbb(j,2))/(2.d0*area(j))       


      uH1=uH1+(w(i)*((gux(i)-guhx(i))**2+(guy(i)-guhy(i))**2))

      p1(i) = nabla1(i)*(2.d0*nabla1(i)-1.d0)
      p2(i) = nabla2(i)*(2.d0*nabla2(i)-1.d0)
      p3(i) = nabla3(i)*(2.d0*nabla3(i)-1.d0)
      p4(i) = 4.d0*nabla2(i)*nabla3(i)
      p5(i) = 4.d0*nabla1(i)*nabla3(i)
      p6(i) = 4.d0*nabla1(i)*nabla2(i)

5100  CONTINUE
      temp_M(1,1)=SUM(w*p1*p1)
      temp_M(1,2)=SUM(w*p1*p2)
      temp_M(1,3)=SUM(w*p1*p3)
      temp_M(1,4)=SUM(w*p1*p4)
      temp_M(1,5)=SUM(w*p1*p5)
      temp_M(1,6)=SUM(w*p1*p6)

      temp_M(2,1)=SUM(w*p2*p1)
      temp_M(2,2)=SUM(w*p2*p2)
      temp_M(2,3)=SUM(w*p2*p3)
      temp_M(2,4)=SUM(w*p2*p4)
      temp_M(2,5)=SUM(w*p2*p5)
      temp_M(2,6)=SUM(w*p2*p6)

      temp_M(3,1)=SUM(w*p3*p1)
      temp_M(3,2)=SUM(w*p3*p2)
      temp_M(3,3)=SUM(w*p3*p3)
      temp_M(3,4)=SUM(w*p3*p4)
      temp_M(3,5)=SUM(w*p3*p5)
      temp_M(3,6)=SUM(w*p3*p6)

      temp_M(4,1)=SUM(w*p4*p1)
      temp_M(4,2)=SUM(w*p4*p2)
      temp_M(4,3)=SUM(w*p4*p3)
      temp_M(4,4)=SUM(w*p4*p4)
      temp_M(4,5)=SUM(w*p4*p5)
      temp_M(4,6)=SUM(w*p4*p6)

      temp_M(5,1)=SUM(w*p5*p1)
      temp_M(5,2)=SUM(w*p5*p2)
      temp_M(5,3)=SUM(w*p5*p3)
      temp_M(5,4)=SUM(w*p5*p4)
      temp_M(5,5)=SUM(w*p5*p5)
      temp_M(5,6)=SUM(w*p5*p6)

      temp_M(6,1)=SUM(w*p6*p1)
      temp_M(6,2)=SUM(w*p6*p2)
      temp_M(6,3)=SUM(w*p6*p3)
      temp_M(6,4)=SUM(w*p6*p4)
      temp_M(6,5)=SUM(w*p6*p5)
      temp_M(6,6)=SUM(w*p6*p6)


      DO 5010 i=1,quad_num
      temp_b(1)=temp_b(1)+w(i)*Exact_u(xx(i),yy(i),testproblem)*p1(i)
      temp_b(2)=temp_b(2)+w(i)*Exact_u(xx(i),yy(i),testproblem)*p2(i)
      temp_b(3)=temp_b(3)+w(i)*Exact_u(xx(i),yy(i),testproblem)*p3(i)
      temp_b(4)=temp_b(4)+w(i)*Exact_u(xx(i),yy(i),testproblem)*p4(i)
      temp_b(5)=temp_b(5)+w(i)*Exact_u(xx(i),yy(i),testproblem)*p5(i)
      temp_b(6)=temp_b(6)+w(i)*Exact_u(xx(i),yy(i),testproblem)*p6(i)
5010  CONTINUE
      CALL inverse(temp_M,invtemp_M,6)
      DO 5050 k=1,6
      DO 5050 l=1,6
      uu(k)=uu(k)+invtemp_M(k,l)*temp_b(l)
5050  CONTINUE
      uuu=uu(1)*p1+uu(2)*p2+uu(3)*p3+uu(4)*p4+uu(5)*p5+uu(6)*p6
      uQL2 = uQL2+SUM(w*(uuu-uh)**2)

5000  CONTINUE

      uL2 = DSQRT(uL2)
      uH1 = DSQRT(uH1)
      uQL2 = DSQRT(uQL2)

      uerror = 0.d0
      
      DO 5200 i=1,NO+3*NE
      DO 5200 j=1,NO+3*NE
      uerror=uerror+(x(i)-ExactU(i))*AA(i,j)*(x(j)-ExactU(j))
5200  CONTINUE
      uerror = DSQRT(uerror)

      allocate(length(NE))
      DO 5300 i=1,NE
      length(i)=length(i)+(node(edge(i,2),2)-node(edge(i,1),2))**2+ &
                    (node(edge(i,2),1)-node(edge(i,1),1))**2
      length(i) = DSQRT(length(i))
5300  CONTINUE

      ubL2=0.d0
      DO 5400 i=1,NE
      k=NO+NE+2*i-1
      ubL2=ubL2+((length(i)**2)*((x(k)-ExactU(k))**2)/2.d0+ &
                (length(i)**2)*((x(k)-ExactU(k))**2)/2.d0)
5400  CONTINUE
      ubL2 = DSQRT(ubL2)
!-----------Plot out-----------------------------
      CALL PLOT(x,ExactU)
      RETURN
      END
!------------------------------------------------

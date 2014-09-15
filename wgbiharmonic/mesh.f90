!------------------------------------------------
      SUBROUTINE MESH

!   For enumeration number of edges in the paper
!   'Three matlab implementations of the lowest-order
!   Raviart-Thomas MFEM with a posteriori error control'
!   by C.Bahriawati and C.Carstensen

      USE SHAREDATA

      IMPLICIT NONE

      integer :: sizeII,p,p1,p2,p3
      integer :: B(maxnumber,maxnumber)   ! for node2ele
      real(8) :: mid(maxnumber,2)
      real(8) :: temp1,temp2,temp3
       

!  To enumerate number of edges

      node2ele = 0.0       ! Initialization
      temp1 = 0.0;temp2=0.0;temp3=0.0;

      DO 2000 j=1,n_elem 
      node2ele(elem(j,1),elem(j,2))=node2ele(elem(j,1),elem(j,2))+j
      node2ele(elem(j,2),elem(j,3))=node2ele(elem(j,2),elem(j,3))+j
      node2ele(elem(j,3),elem(j,1))=node2ele(elem(j,3),elem(j,1))+j
2000  CONTINUE


      DO 2020 i=1,n_elem
      DO 2020 j=1,n_elem

      B(i,j) = node2ele(i,j)+node2ele(j,i)
2020  CONTINUE

!   Find the non-zero number in upper triangle B matrix
      
      k=0
      DO 2040 j=1,n_elem
      DO 2050 i=1,n_elem

      IF (i.LE.j) THEN
      IF (B(i,j).NE.0) THEN
      k = k+1
      II(k) = i
      JJ(k) = j

      ENDIF
      ENDIF

2050  CONTINUE
2040  CONTINUE

      sizeII = k  ! size of II

      node2edge = 0  ! initialization of node2edge

      DO 2060 i=1,sizeII

      node2edge(II(i),JJ(i)) = i

2060  CONTINUE

      DO 2070 i=1,sizeII
      DO 2080 j=1,sizeII
      
      node2edge(i,j) = node2edge(i,j)+node2edge(j,i)


2080  CONTINUE
2070  CONTINUE
      n_edge = sizeII   ! number of edges

      edge2ele = 0

      DO 2100 i=1,n_elem
      DO 2110 j=1,3

      p = node2edge(elem(i,j),elem(i,MOD(j,3)+1))
      IF(edge2ele(p,1).EQ.0) THEN
        edge2ele(p,1) = elem(i,j)
        edge2ele(p,2) = elem(i,MOD(j,3)+1)
        edge2ele(p,3) = node2ele(elem(i,j),elem(i,MOD(j,3)+1)) 
        edge2ele(p,4) = node2ele(elem(i,MOD(j,3)+1),elem(i,j))
      ENDIF
2110  CONTINUE
2100  CONTINUE

!  interioredge  searching

      k=0
      DO i=1,n_edge
      IF(edge2ele(i,4).NE.0) THEN
      k=k+1
      DO j=1,4
      interioredge(k,j) = edge2ele(i,j)
      ENDDO
      ENDIF
      ENDDO

      n_in = k

!  exterioredges
      II = 0
      JJ = 0
      k=0
      DO i=1,n_edge
      IF(edge2ele(i,4)==0)  THEN
      k=k+1
      II(k) = i
      JJ(k) = i
      ENDIF
      ENDDO

      k=0
      DO i=1,n_edge
      
      IF(edge2ele(i,4)==0) THEN
      k=k+1
      DO j=1,3
      exterioredge(k,j) = edge2ele(i,j)
      ENDDO
      exterioredge(k,4) = II(k)
      ENDIF
      ENDDO

      n_ex = k

!  Connectivity of element to edge
      ele2edge = 0.0
      
      DO 2200 i=1,n_elem
      p1 = elem(i,1)
      p2 = elem(i,2)
      p3 = elem(i,3)
      
      ele2edge(i,1) = node2edge(p2,p3)
      ele2edge(i,2) = node2edge(p1,p3)
      ele2edge(i,3) = node2edge(p2,p1)
2200  CONTINUE

      edge = 0.0   ! edge matrix initialzation
      
      DO 2300 i=1,n_edge
      DO 2300 j=1,2
      edge(i,j) = edge2ele(i,j)
2300  CONTINUE
    
      n_mid = n_edge

      mid = 0.0d0   ! mid matrix initialization

      DO 2400 i=1,n_mid
      mid(i,1) = (node(edge(i,1),1)+node(edge(i,2),1))/2.0d0
      mid(i,2) = (node(edge(i,1),2)+node(edge(i,2),2))/2.0d0
2400  CONTINUE

      ele2mid = ele2edge

      RETURN 
      END
!  End of connectivity
!------------------------------------------------

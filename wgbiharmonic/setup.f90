
!------------------------------------------------
      MODULE SHAREDATA
      
      IMPLICIT NONE
      SAVE

      integer :: maxnumber
      parameter (maxnumber = 5000)
      integer :: testproblem  ! parameter
      integer :: i,j,k,l,m
      integer :: iter,n,nn,n_elem,n_node,n_edge,n_mid
      real(8) :: node(maxnumber,2)
      integer :: elem(maxnumber,3)
!
      integer :: node2ele(maxnumber,maxnumber)
      integer :: edge2ele(maxnumber,4)
      integer :: II(maxnumber),JJ(maxnumber) ! position of non-zero value in matrix B
      integer :: node2edge(maxnumber,maxnumber)
      integer :: ele2edge(maxnumber,3),edge(maxnumber,2)
      integer :: ele2mid(maxnumber,3)
      integer :: exterioredge(maxnumber,4),interioredge(maxnumber,4)
      integer :: n_ex,n_in;
      integer :: quad_num

      real(8) :: uL2,uH1,uQL2,uerror,ubL2,u0Inf
      
      real(8) :: QL2error(100),L2inf(100),H1error(100)
      real(8) :: Newerror(100),L2error(100),L2berror(100)
      real(8) :: Con(100),hh(100)
      

      END MODULE SHAREDATA
!------------------------------------------------

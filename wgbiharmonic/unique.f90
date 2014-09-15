!------------------------------------------------
      SUBROUTINE Unique(arraya,numrows,BNode)

      USE SHAREDATA
      IMPLICIT NONE

      INTEGER :: numrows

      INTEGER, DIMENSION(numrows) :: arraya
      INTEGER, DIMENSION(numrows) :: BNode
      LOGICAL, DIMENSION(:), ALLOCATABLE :: mask
      INTEGER :: ix
      INTEGER, DIMENSION(:), ALLOCATABLE :: index_vector


        ! First, find the duplicate elements
      ALLOCATE(mask(numrows))
      mask = .TRUE.

      DO ix = numrows,2,-1
      mask(ix) = .NOT.(ANY(arraya(:ix-1)==arraya(ix)))
      END DO

      ! Make an index vector
      ALLOCATE(index_vector, source=PACK([(ix, ix=1,numrows) ],mask))

     ! Now copy the unique elements of a into b
      BNode=arraya(index_vector)
      CALL Sort(BNode,size(BNode))

      RETURN
      END
!------------------------------------------------
      INTEGER FUNCTION FindMinimum(x,Start,End)
!  This function returns the location of the minimum in the 
!  section between Start and End
       
      IMPLICIT NONE
      INTEGER, INTENT(IN)                :: Start,End
      INTEGER, DIMENSION(End)            :: x
      INTEGER                            :: Minimum
      INTEGER                            :: location
      INTEGER                            :: i

      Minimum = x(Start)          ! asume the first is the min
      location = Start            ! record its position 
      DO i=Start+1,End            ! start with next elements
      IF (x(i)<Minimum) THEN      ! if x(i) less than the min
      Minimum = x(i)              ! Yes, a new minimum found
      location = i                ! record its position
      END IF
      END DO
      FindMinimum = location     ! return the position
      END FUNCTION FindMinimum 
!------------------------------------------------
      SUBROUTINE Sort(x,size)
!  This subroutine receive an array x() and sorts
!  it into ascending order

      IMPLICIT NONE
      INTEGER                               :: size
      INTEGER, DIMENSION(size)              :: x
      INTEGER                               :: location
      INTEGER                               :: i
      INTEGER                               :: FindMinimum

      DO i=1,size-1            ! except for the last
      location = FindMinimum(x,i,size)   ! find min from this to last
      CALL Swap(x(i),x(location))  ! swap this and the minimum
      END DO
      END SUBROUTINE Sort
!------------------------------------------------
      SUBROUTINE Swap(a,b)
!  This subroutine swaps the values of its two formal 
!  arguments

      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: a,b
      INTEGER               :: Temp

      Temp = a
      a = b
      b = Temp
      END SUBROUTINE Swap
!------------------------------------------------

!------------------------------------------------
       SUBROUTINE QUAD_RULE(quad_number,quad_xy,quad_w)
! --- focus on quad_num=13 

       IMPLICIT NONE

       REAL(8) :: a,b,c,d,e,f,g,h,w,t,u,v
       integer :: quad_number,j
       REAL(8) :: quad_xy(quad_number,2),quad_w(quad_number)

       IF (quad_number.EQ.1) THEN
!       quad_xy(1,1) = 0.5773502691896257645091488
!       quad_xy(1,2) = 0.5773502691896257645091488
       quad_w(1) = 1.d0

       ELSE IF(quad_number.EQ.3) THEN

       DO j=1,2
!       quad_xy(1,j) = 0.d0
!       quad_xy(2,j) = 0.77459666
!       quad_xy(3,j) = 0.77459666
       ENDDO
       quad_w = 1.d0/3.d0

       ELSE IF(quad_number.EQ.13) THEN

       h = 1.d0/3.d0
       a = 0.479308067841923
       b = 0.260345966079038
       c = 0.869739794195568
       d = 0.065130102902216
       e = 0.638444188569809
       f = 0.312865496004875
       g = 0.048690315425316
        
       w = -0.149570044467670
       t =  0.175615257433204
       u =  0.053347235608839
       v =  0.077113760890257

       quad_xy(1,1) = h
       quad_xy(1,2) = h
       quad_xy(2,1) = a
       quad_xy(2,2) = b
       quad_xy(3,1) = b
       quad_xy(3,2) = a
       quad_xy(4,1) = b
       quad_xy(4,2) = b
       quad_xy(5,1) = c
       quad_xy(5,2) = d
       quad_xy(6,1) = d
       quad_xy(6,2) = c
       quad_xy(7,1) = d
       quad_xy(7,2) = d
       quad_xy(8,1) = e
       quad_xy(8,2) = f
       quad_xy(9,1) = e
       quad_xy(9,2) = g
       quad_xy(10,1) = f
       quad_xy(10,2) = e
       quad_xy(11,1) = f
       quad_xy(11,2) = g
       quad_xy(12,1) = g
       quad_xy(12,2) = e
       quad_xy(13,1) = g
       quad_xy(13,2) = f

       quad_w(1) = w
       quad_w(2) = t
       quad_w(3) = t
       quad_w(4) = t
       quad_w(5) = u
       quad_w(6) = u
       quad_w(7) = u
       quad_w(8) = v
       quad_w(9) = v
       quad_w(10) = v
       quad_w(11) = v
       quad_w(12) = v
       quad_w(13) = v
       ENDIF

       RETURN
       END
!------------------------------------------------

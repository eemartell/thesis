#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)    
    ! Declarations
	implicit none

    ! mexFunction argument
    integer*4, intent(in) :: nlhs, nrhs
    mwpointer, intent(in), dimension(*) :: prhs
    mwpointer, intent(out), dimension(*) :: plhs
    
    ! Function declarations
    integer, pointer :: mxGetPr
    mwsize :: mxCreateDoubleMatrix,mxGetM,mxGetN
    double precision :: mxGetScalar
    
    ! Inputs
	mwsize nx1,nx2,nx3,nx4,nx5,nx6,nx7,nodes
	integer, pointer :: x1,x2,x3,x4,x5,x6,x7
    integer, pointer :: x1i,x2i,x3i,x4i,x5i,x6i,x7i
	integer, pointer :: r1,r2,r3,r4
	
	! Outputs
	integer, pointer :: o1,o2,o3,o4
           
    ! Load Inputs
    ! Number of Points
    nx1 = mxGetN(prhs(1))
    nx2 = mxGetN(prhs(2))
    nx3 = mxGetN(prhs(3))
    nx4 = mxGetN(prhs(4))
    nx5 = mxGetN(prhs(5))
    nx6 = mxGetN(prhs(6))
    nx7 = mxGetN(prhs(7))
    
    ! Grid Values
    x1 => mxGetPr(prhs(1))
    x2 => mxGetPr(prhs(2))
    x3 => mxGetPr(prhs(3))
    x4 => mxGetPr(prhs(4))
    x5 => mxGetPr(prhs(5))
    x6 => mxGetPr(prhs(6))
    x7 => mxGetPr(prhs(7))
    
    ! Point to evaluate
    x1i => mxGetPr(prhs(8))
    x2i => mxGetPr(prhs(9))
    x3i => mxGetPr(prhs(10))
    x4i => mxGetPr(prhs(11))
    nodes = mxGetM(prhs(12))
    x5i => mxGetPr(prhs(12))
    x6i => mxGetPr(prhs(13))
    x7i => mxGetPr(prhs(14))
    
    ! Rules
    r1 => mxGetPr(prhs(15))
    r2 => mxGetPr(prhs(16))
    r3 => mxGetPr(prhs(17))
    r4 => mxGetPr(prhs(18))
    
    !Create array for return argument
    plhs(1) = mxCreateDoubleMatrix(nodes,1,0)   
    plhs(2) = mxCreateDoubleMatrix(nodes,1,0)   
    plhs(3) = mxCreateDoubleMatrix(nodes,1,0)
    plhs(4) = mxCreateDoubleMatrix(nodes,1,0)
    o1 => mxGetPr(plhs(1))
    o2 => mxGetPr(plhs(2))
    o3 => mxGetPr(plhs(3))
    o4 => mxGetPr(plhs(4))
    
    ! Call subroutine for assignment
    call allterp(	o1,o2,o3,o4, &
					nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
					x1,x2,x3,x4,x5,x6,x7, &
					x1i,x2i,x3i,x4i,x5i,x6i,x7i, &
					r1,r2,r3,r4,nodes)

end subroutine mexFunction

subroutine allterp(	o1,o2,o3,o4, &
					nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
					x1,x2,x3,x4,x5,x6,x7, &
					x1i,x2i,x3i,x4i,x5i,x6i,x7i, &
					r1,r2,r3,r4,nodes)

    implicit none
    mwsize nx1,nx2,nx3,nx4,nx5,nx6,nx7,nodes
    mwsize i
    double precision fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4
    double precision, dimension(nx5) :: x5
    double precision, dimension(nx6) :: x6
    double precision, dimension(nx7) :: x7
    double precision, dimension(nx1,nx2,nx3,nx4,nx5,nx6,nx7) :: r1,r2,r3,r4
    double precision, dimension(nodes) :: x1i,x2i,x3i,x4i,x5i,x6i,x7i
    double precision, dimension(nodes) :: o1,o2,o3,o4
	  
    ! Interpolate
    do i = 1,nodes  
        o1(i) = fastrap(	nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
							x1,x2,x3,x4,x5,x6,x7, &
							x1i(i),x2i(i),x3i(i),x4i(i),x5i(i),x6i(i),x7i(i), &
							r1)
        o2(i) = fastrap(	nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
							x1,x2,x3,x4,x5,x6,x7, &
							x1i(i),x2i(i),x3i(i),x4i(i),x5i(i),x6i(i),x7i(i), &
							r2)
        o3(i) = fastrap(	nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
							x1,x2,x3,x4,x5,x6,x7, &
							x1i(i),x2i(i),x3i(i),x4i(i),x5i(i),x6i(i),x7i(i), &
							r3)
        o4(i) = fastrap(	nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
							x1,x2,x3,x4,x5,x6,x7, &
							x1i(i),x2i(i),x3i(i),x4i(i),x5i(i),x6i(i),x7i(i), &
							r4)
    end do
end subroutine allterp


function fastrap(	nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
					x1,x2,x3,x4,x5,x6,x7, &
					x1i,x2i,x3i,x4i,x5i,x6i,x7i,z)

    implicit none
    mwSize nx1,nx2,nx3,nx4,nx5,nx6,nx7
    double precision :: x1i,x2i,x3i,x4i,x5i,x6i,x7i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4
    double precision, dimension(nx5) :: x5
    double precision, dimension(nx6) :: x6
    double precision, dimension(nx7) :: x7
    double precision, dimension(nx1,nx2,nx3,nx4,nx5,nx6,nx7) :: z
    
    double precision :: s1, s2, s3, s4, s5, s6,s7
    double precision :: x1i_min, x2i_min, x3i_min, x4i_min, x5i_min, x6i_min, x7i_min
    mwSize loc1, loc2, loc3, loc4, loc5, loc6, loc7
    double precision, dimension(7) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2, w3, w4, w5, w6, w7
    mwSize m1, m2, m3, m4, m5, m6, m7
    
    s1 = x1(2) - x1(1)
    x1i_min = x1i - x1(1)
    loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));
    
    s2 = x2(2) - x2(1)
    x2i_min = x2i - x2(1)
    loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));
    
    s3 = x3(2) - x3(1)
    x3i_min = x3i - x3(1)
    loc3 = min(nx3-1,max(1,floor(x3i_min/s3) + 1));
    
    s4 = x4(2) - x4(1)
    x4i_min = x4i - x4(1)
    loc4 = min(nx4-1,max(1,floor(x4i_min/s4) + 1));
    
    s5 = x5(2) - x5(1)
    x5i_min = x5i - x5(1)
    loc5 = min(nx5-1,max(1,floor(x5i_min/s5) + 1));   
    
    s6 = x6(2) - x6(1)
    x6i_min = x6i - x6(1)
    loc6 = min(nx6-1,max(1,floor(x6i_min/s6) + 1));  
    
    s7 = x7(2) - x7(1)
    x7i_min = x7i - x7(1)
    loc7 = min(nx7-1,max(1,floor(x7i_min/s7) + 1));             

    xi = [x1i, x2i, x3i, x4i, x5i, x6i, x7i]
    xi_left = [x1(loc1), x2(loc2), x3(loc3), x4(loc4), x5(loc5), x6(loc6), x7(loc7)]
    xi_right = [x1(loc1+1), x2(loc2+1), x3(loc3+1), x4(loc4+1), x5(loc5+1), x6(loc6+1), x7(loc7+1)]

    w_2 = (xi - xi_left)/(xi_right - xi_left)
    w_1 = 1 - w_2
    w1 = [w_1(1), w_2(1)]
    w2 = [w_1(2), w_2(2)]
    w3 = [w_1(3), w_2(3)]
    w4 = [w_1(4), w_2(4)]
    w5 = [w_1(5), w_2(5)]
    w6 = [w_1(6), w_2(6)]
    w7 = [w_1(7), w_2(7)]
    
    fastrap = 0
    
    do m7 = 0, 1
     do m6 = 0, 1
      do m5 = 0, 1
       do m4 = 0, 1
        do m3 = 0, 1
         do m2 = 0, 1
          do m1 = 0, 1  
           fastrap = fastrap + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*w6(m6+1)*w7(m7+1)* &
							z(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5,loc6+m6,loc7+m7)
          end do
		 end do
	    end do  
       end do
      end do
     end do
    end do

end function fastrap
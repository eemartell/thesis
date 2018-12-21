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
	mwsize nx1,nx2,nx3,nodes
	integer, pointer :: x1,x2,x3
    integer, pointer :: x1i,x2i,x3i
	integer, pointer :: r1,r2
	
	! Outputs
	integer, pointer :: o1,o2
           
    ! Load Inputs
    ! Number of Points
    nx1 = mxGetN(prhs(1))
    nx2 = mxGetN(prhs(2))
    nx3 = mxGetN(prhs(3))
    
    ! Grid Values
    x1 => mxGetPr(prhs(1))
    x2 => mxGetPr(prhs(2))
    x3 => mxGetPr(prhs(3))
    
    ! Point to evaluate
    nodes = mxGetM(prhs(4))
    x1i => mxGetPr(prhs(4))
    x2i => mxGetPr(prhs(5))
    x3i => mxGetPr(prhs(6))
    
    ! Rules
    r1 => mxGetPr(prhs(7))
    r2 => mxGetPr(prhs(8))
    
    !Create array for return argument
    plhs(1) = mxCreateDoubleMatrix(nodes,1,0)   
    plhs(2) = mxCreateDoubleMatrix(nodes,1,0)   
    o1 => mxGetPr(plhs(1))   
    o2 => mxGetPr(plhs(2))   
    
    ! Call subroutine for assignment
    call allterp(o1,o2,nx1,nx2,nx3,x1,x2,x3,x1i,x2i,x3i,r1,r2,nodes)

end subroutine mexFunction

subroutine allterp(o1,o2,nx1,nx2,nx3,x1,x2,x3,x1i,x2i,x3i,r1,r2,nodes)

    implicit none
    mwsize nx1,nx2,nx3,nodes
    mwsize i
    double precision fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx1,nx2,nx3) :: r1,r2
    double precision, dimension(nodes) :: x1i,x2i,x3i
    double precision, dimension(nodes) :: o1,o2

    ! Interpolate
    do i = 1,nodes  
        o1(i) = fastrap(nx1,nx2,nx3,x1,x2,x3,x1i(i),x2i(i),x3i(i),r1)
        o2(i) = fastrap(nx1,nx2,nx3,x1,x2,x3,x1i(i),x2i(i),x3i(i),r2)
    end do
end subroutine allterp


function fastrap(nx1,nx2,nx3,x1,x2,x3,x1i,x2i,x3i,z)

    implicit none
    mwsize nx1,nx2,nx3
    double precision :: x1i,x2i,x3i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx1,nx2,nx3) :: z
    
    double precision :: s1, s2, s3
    double precision :: x1i_min, x2i_min, x3i_min
    mwsize loc1, loc2, loc3
    double precision, dimension(3) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2, w3
    mwsize m1, m2, m3
    
    s1 = x1(2) - x1(1)
    x1i_min = x1i - x1(1)
    loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));
    
    s2 = x2(2) - x2(1)
    x2i_min = x2i - x2(1)
    loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));
    
    s3 = x3(2) - x3(1)
    x3i_min = x3i - x3(1)
    loc3 = min(nx3-1,max(1,floor(x3i_min/s3) + 1));     
    
    xi = [x1i, x2i, x3i]
    xi_left = [x1(loc1), x2(loc2), x3(loc3)]
    xi_right = [x1(loc1+1), x2(loc2+1), x3(loc3+1)]

    w_2 = (xi - xi_left)/(xi_right - xi_left)
    w_1 = 1 - w_2
    w1 = [w_1(1), w_2(1)]
    w2 = [w_1(2), w_2(2)]
    w3 = [w_1(3), w_2(3)]
    
    fastrap = 0
    
    do m3 = 0, 1
     do m2 = 0, 1
      do m1 = 0, 1
       fastrap = fastrap + w1(m1+1)*w2(m2+1)*w3(m3+1)*z(loc1+m1,loc2+m2,loc3+m3)
      end do
     end do
    end do
    
end function fastrap
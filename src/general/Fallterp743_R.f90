#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    ! Declarations
	implicit none

    ! mexFunction argument
    mwPointer plhs(*), prhs(*)    
    integer*4 nlhs, nrhs
    
    ! Function declarations
    mwSize mxGetN  
    mwpointer mxGetPr, mxCreateNumericArray, mxGetDimensions 
    double precision mxGetScalar  
    integer*4 mxClassIDFromClassName 
    
    ! Pointers to input/output mxArrays
    mwpointer x4_pr,x5_pr,x6_pr,x7_pr
    mwpointer pf1_pr,pf2_pr,pf3_pr,pf4_pr
    mwpointer o1_pr,o2_pr,o3_pr,o4_pr
    
    ! Array information
	mwSize nx1,nx2,nx3,nx4,nx5,nx6,nx7,nodes
	integer*4 myclassid
	double precision, allocatable, dimension(:) :: x4,x5,x6,x7
	double precision x4i,x5i,x6i,x7i
	double precision, allocatable, dimension(:,:,:,:,:,:,:) :: pf1,pf2,pf3,pf4
	double precision, allocatable, dimension(:,:,:) :: o1,o2,o3,o4
	  
    ! Load Inputs
    nx1 = mxGetScalar(prhs(1))
    nx2 = mxGetScalar(prhs(2))
    nx3 = mxGetScalar(prhs(3))
    nx4 = mxGetScalar(prhs(4))
    nx5 = mxGetScalar(prhs(5))
    nx6 = mxGetScalar(prhs(6))
    nx7 = mxGetScalar(prhs(7))
    nodes = nx1*nx2*nx3*nx4*nx5*nx6*nx7
    allocate(x4(nx4))
    allocate(x5(nx5))
    allocate(x6(nx6))
    allocate(x7(nx7))
    x4_pr = mxGetPr(prhs(8))
    x5_pr = mxGetPr(prhs(9))
    x6_pr = mxGetPr(prhs(10))
    x7_pr = mxGetPr(prhs(11))
    call mxCopyPtrToReal8(x4_pr,x4,nx4)
    call mxCopyPtrToReal8(x5_pr,x5,nx5)
    call mxCopyPtrToReal8(x6_pr,x6,nx6)
    call mxCopyPtrToReal8(x7_pr,x7,nx7)
    ! Point to evaluate
    x4i = mxGetScalar(prhs(12))
    x5i = mxGetScalar(prhs(13))
    x6i = mxGetScalar(prhs(14))
    x7i = mxGetScalar(prhs(15))
	  
    ! Policy functions
    allocate(pf1(nx1,nx2,nx3,nx4,nx5,nx6,nx7))
    allocate(pf2(nx1,nx2,nx3,nx4,nx5,nx6,nx7))
    allocate(pf3(nx1,nx2,nx3,nx4,nx5,nx6,nx7))
    allocate(pf4(nx1,nx2,nx3,nx4,nx5,nx6,nx7))
    pf1_pr = mxGetPr(prhs(16))
    pf2_pr = mxGetPr(prhs(17))
    pf3_pr = mxGetPr(prhs(18))
    pf4_pr = mxGetPr(prhs(19))
    call mxCopyPtrToReal8(pf1_pr,pf1,nodes)
    call mxCopyPtrToReal8(pf2_pr,pf2,nodes)
    call mxCopyPtrToReal8(pf3_pr,pf3,nodes)
    call mxCopyPtrToReal8(pf4_pr,pf4,nodes)
	  
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(o1(nx1,nx2,nx3))
    allocate(o2(nx1,nx2,nx3))
    allocate(o3(nx1,nx2,nx3))
    allocate(o4(nx1,nx2,nx3))
    plhs(1) = mxCreateNumericArray(3,[nx1,nx2,nx3],myclassid,0)
    plhs(2) = mxCreateNumericArray(3,[nx1,nx2,nx3],myclassid,0)
    plhs(3) = mxCreateNumericArray(3,[nx1,nx2,nx3],myclassid,0)
    plhs(4) = mxCreateNumericArray(3,[nx1,nx2,nx3],myclassid,0)
    o1_pr = mxGetPr(plhs(1))
    o2_pr = mxGetPr(plhs(2))
    o3_pr = mxGetPr(plhs(3))
    o4_pr = mxGetPr(plhs(4))
                      
    ! Call subroutine for assignment
    call allterp(	nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
					x4,x5,x6,x7, &
					x4i,x5i,x6i,x7i, &
					pf1,pf2,pf3,pf4, &
					o1,o2,o3,o4)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(o1,o1_pr,nx1*nx2*nx3)
    call mxCopyReal8ToPtr(o2,o2_pr,nx1*nx2*nx3)
    call mxCopyReal8ToPtr(o3,o3_pr,nx1*nx2*nx3)
    call mxCopyReal8ToPtr(o4,o4_pr,nx1*nx2*nx3)
    
    ! Deallocate arrays
    deallocate(x4)
    deallocate(x5)
    deallocate(x6)
    deallocate(x7)
    deallocate(pf1)
    deallocate(pf2)
    deallocate(pf3)
    deallocate(pf4)
    deallocate(o1) 
    deallocate(o2) 
    deallocate(o3) 
    deallocate(o4)
    
end subroutine mexFunction

subroutine allterp(	nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
					x4,x5,x6,x7, &
					x4i,x5i,x6i,x7i, &
					pf1,pf2,pf3,pf4, &
					o1,o2,o3,o4)

    implicit none
    mwSize :: nx1,nx2,nx3,nx4,nx5,nx6,nx7
    double precision :: x4i,x5i,x6i,x7i,x8i,x4(nx4),x5(nx5),x6(nx6),x7(nx7)
    double precision, dimension(nx1,nx2,nx3,nx4,nx5,nx6,nx7) :: pf1,pf2,pf3,pf4
    double precision, dimension(nx1,nx2,nx3) :: o1,o2,o3,o4
    
    double precision :: s4, s5, s6, s7
    double precision :: x4i_min, x5i_min, x6i_min, x7i_min
    mwSize loc4, loc5, loc6, loc7
    double precision, dimension(4) :: xi, xi_left, xi_right, w_2, w_1
	double precision, dimension(2) :: w4, w5, w6, w7
	mwSize :: m4, m5, m6, m7
	double precision :: wtemp

	s4 = x4(2) - x4(1)
	s5 = x5(2) - x5(1)
	s6 = x6(2) - x6(1)
	s7 = x7(2) - x7(1)
	
	x4i_min = x4i - x4(1)
	loc4 = min(nx4-1,max(1,floor(x4i_min/s4) + 1)); 
	
	x5i_min = x5i - x5(1)
	loc5 = min(nx5-1,max(1,floor(x5i_min/s5) + 1)); 
	
	x6i_min = x6i - x6(1)
	loc6 = min(nx6-1,max(1,floor(x6i_min/s6) + 1));
	
	x7i_min = x7i - x7(1)
	loc7 = min(nx7-1,max(1,floor(x7i_min/s7) + 1));

	xi = [x4i, x5i, x6i, x7i]
	xi_left = [x4(loc4), x5(loc5), x6(loc6), x7(loc7)]
	xi_right = [x4(loc4+1), x5(loc5+1), x6(loc6+1), x7(loc7+1)]

	w_2 = (xi - xi_left)/(xi_right - xi_left)
	w_1 = 1 - w_2
	
	w4 = [w_1(1), w_2(1)]
	w5 = [w_1(2), w_2(2)]
	w6 = [w_1(3), w_2(3)]
	w7 = [w_1(4), w_2(4)]
	
	o1 = 0.d0
	o2 = 0.d0
	o3 = 0.d0
	o4 = 0.d0
	  do m7 = 0,1
		do m6 = 0,1
		  do m5 = 0,1
		   do m4 = 0,1
			wtemp = w4(m4+1)*w5(m5+1)*w6(m6+1)*w7(m7+1)
			o1 = o1 + wtemp*pf1(:,:,:,loc4+m4,loc5+m5,loc6+m6,loc7+m7)
			o2 = o2 + wtemp*pf2(:,:,:,loc4+m4,loc5+m5,loc6+m6,loc7+m7)
			o3 = o3 + wtemp*pf3(:,:,:,loc4+m4,loc5+m5,loc6+m6,loc7+m7)
			o4 = o4 + wtemp*pf4(:,:,:,loc4+m4,loc5+m5,loc6+m6,loc7+m7)
		   end do
		  end do
		end do
	  end do
			  
end subroutine allterp
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
    mwpointer x4_pr
    mwpointer pf1_pr,pf2_pr
    mwpointer o1_pr,o2_pr
    
    ! Array information
	mwSize nx1,nx2,nx3,nx4,nx5,nodes
	integer*4 myclassid
	double precision, allocatable, dimension(:) :: x4
	double precision x4i
	double precision, allocatable, dimension(:,:,:,:) :: pf1,pf2
	double precision, allocatable, dimension(:,:,:) :: o1,o2
	  
    ! Load Inputs
    nx1 = mxGetScalar(prhs(1))
    nx2 = mxGetScalar(prhs(2))
    nx3 = mxGetScalar(prhs(3))
    nx4 = mxGetScalar(prhs(4))
    nodes = nx1*nx2*nx3*nx4
    allocate(x4(nx4))
    x4_pr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(x4_pr,x4,nx4) 
    ! Point to evaluate
    x4i = mxGetScalar(prhs(6))
	  
    ! Policy functions
    allocate(pf1(nx1,nx2,nx3,nx4))
    allocate(pf2(nx1,nx2,nx3,nx4))
    pf1_pr = mxGetPr(prhs(7))
    pf2_pr = mxGetPr(prhs(8))
    call mxCopyPtrToReal8(pf1_pr,pf1,nodes)
    call mxCopyPtrToReal8(pf2_pr,pf2,nodes)
	  
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(o1(nx1,nx2,nx3)) 
    allocate(o2(nx1,nx2,nx3))
    plhs(1) = mxCreateNumericArray(3,[nx1,nx2,nx3],myclassid,0)
    plhs(2) = mxCreateNumericArray(3,[nx1,nx2,nx3],myclassid,0)
    o1_pr = mxGetPr(plhs(1))      
    o2_pr = mxGetPr(plhs(2))
                      
    ! Call subroutine for assignment
    call allterp(nx1,nx2,nx3,nx4,x4,x4i,pf1,pf2,o1,o2)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(o1,o1_pr,nx1*nx2*nx3)
    call mxCopyReal8ToPtr(o2,o2_pr,nx1*nx2*nx3)
    
    ! Deallocate arrays
    deallocate(x4)
    deallocate(pf1)
    deallocate(pf2)    
    deallocate(o1) 
    deallocate(o2)     
    
end subroutine mexFunction

subroutine allterp(	nx1,nx2,nx3,nx4, &
					x4, &
					x4i, &
					pf1,pf2, &
					o1,o2)

    implicit none
    mwSize :: nx1,nx2,nx3,nx4,nx5
    double precision :: x4i,x5i,x4(nx4)
    double precision, dimension(nx1,nx2,nx3,nx4) :: pf1,pf2
    double precision, dimension(nx1,nx2,nx3) :: o1,o2
    
    double precision :: s4
    double precision :: x4i_min
    mwSize loc4
    double precision, dimension(1) :: xi, xi_left, xi_right, w_2, w_1
    double precision :: w11,w12,w21,w22

    s4 = x4(2) - x4(1)
    
    x4i_min = x4i - x4(1)
    loc4 = min(nx4-1,max(1,floor(x4i_min/s4) + 1));   

	xi = [x4i]
	xi_left = [x4(loc4)]
	xi_right = [x4(loc4+1)]

	w_2 = (xi - xi_left)/(xi_right - xi_left)
	w_1 = 1 - w_2
    
    w11 = w_1(1)
    w21 = w_2(1)

	o1 =  	w11*pf1(:,:,:,loc4) &
		  + w21*pf1(:,:,:,loc4+1)
	o2 =  	w11*pf2(:,:,:,loc4) &
		  + w21*pf2(:,:,:,loc4+1)
			  
end subroutine allterp

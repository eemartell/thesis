#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    ! Declarations
	implicit none

    ! mexFunction argument
    mwPointer plhs(*), prhs(*)    
    integer*4 nlhs, nrhs
    
    ! Function declarations
    mwSize mxGetM,mxGetN  
    mwpointer mxGetPr, mxCreateNumericArray, mxGetDimensions 
    double precision mxGetScalar  
    integer*4 mxClassIDFromClassName 
    
    ! Pointers to input/output mxArrays
    mwpointer x1_pr,x2_pr,x3_pr,x4_pr,x5_pr,x6_pr,x7_pr
    mwpointer x5i_pr,x6i_pr,x7i_pr
    mwpointer pf1_pr,pf2_pr,pf3_pr,pf4_pr
    mwpointer o1_pr,o2_pr,o3_pr,o4_pr
    
    ! Array information
	mwSize nx1,nx2,nx3,nx4,nx5,nx6,nx7,nodes,e5,e6,e7
	integer*4 myclassid
	double precision, allocatable, dimension(:) :: x1,x2,x3,x4,x5,x6,x7,x5i,x6i,x7i
	double precision x1i,x2i,x3i,x4i
	double precision, allocatable, dimension(:,:,:,:,:,:,:) :: pf1,pf2,pf3,pf4
	double precision, allocatable, dimension(:,:,:) :: o1,o2,o3,o4
    
    ! Load Inputs
    ! Grids
    nx1 = mxGetN(prhs(1))  
    nx2 = mxGetN(prhs(2))  
    nx3 = mxGetN(prhs(3))
    nx4 = mxGetN(prhs(4))
    nx5 = mxGetN(prhs(5))
    nx6 = mxGetN(prhs(6))
    nx7 = mxGetN(prhs(7))
    nodes = nx1*nx2*nx3*nx4*nx5*nx6*nx7
    allocate(x1(nx1))
    allocate(x2(nx2))
    allocate(x3(nx3))
    allocate(x4(nx4))
    allocate(x5(nx5))
    allocate(x6(nx6))
    allocate(x7(nx7))
    x1_pr = mxGetPr(prhs(1))
    x2_pr = mxGetPr(prhs(2))
    x3_pr = mxGetPr(prhs(3))
    x4_pr = mxGetPr(prhs(4)) 
    x5_pr = mxGetPr(prhs(5)) 
    x6_pr = mxGetPr(prhs(6))
    x7_pr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(x1_pr,x1,nx1)
    call mxCopyPtrToReal8(x2_pr,x2,nx2)
    call mxCopyPtrToReal8(x3_pr,x3,nx3)
    call mxCopyPtrToReal8(x4_pr,x4,nx4)
    call mxCopyPtrToReal8(x5_pr,x5,nx5)
    call mxCopyPtrToReal8(x6_pr,x6,nx6)
    call mxCopyPtrToReal8(x7_pr,x7,nx7)  
    ! Point to evaluate
    x1i = mxGetScalar(prhs(8))  
    x2i = mxGetScalar(prhs(9))  
    x3i = mxGetScalar(prhs(10)) 
    x4i = mxGetScalar(prhs(11))   
    e5 = mxGetM(prhs(12))
    e6 = mxGetM(prhs(13))   
    e7 = mxGetM(prhs(14)) 
    allocate(x5i(e5))
    allocate(x6i(e6))
    allocate(x7i(e7))
    x5i_pr = mxGetPr(prhs(12))
    x6i_pr = mxGetPr(prhs(13))
    x7i_pr = mxGetPr(prhs(14))
    call mxCopyPtrToReal8(x5i_pr,x5i,e5)
    call mxCopyPtrToReal8(x6i_pr,x6i,e6)
    call mxCopyPtrToReal8(x7i_pr,x7i,e7)

    ! Rules
    allocate(pf1(nx1,nx2,nx3,nx4,nx5,nx6,nx7))
    allocate(pf2(nx1,nx2,nx3,nx4,nx5,nx6,nx7))
    allocate(pf3(nx1,nx2,nx3,nx4,nx5,nx6,nx7))
    allocate(pf4(nx1,nx2,nx3,nx4,nx5,nx6,nx7))
    pf1_pr = mxGetPr(prhs(15))
    pf2_pr = mxGetPr(prhs(16))
    pf3_pr = mxGetPr(prhs(17))
    pf4_pr = mxGetPr(prhs(18))
    call mxCopyPtrToReal8(pf1_pr,pf1,nodes)
    call mxCopyPtrToReal8(pf2_pr,pf2,nodes)
    call mxCopyPtrToReal8(pf3_pr,pf3,nodes)
    call mxCopyPtrToReal8(pf4_pr,pf4,nodes)

    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(o1(e5,e6,e7)) 
    allocate(o2(e5,e6,e7)) 
    allocate(o3(e5,e6,e7)) 
    allocate(o4(e5,e6,e7))
    plhs(1) = mxCreateNumericArray(3,[e5,e6,e7],myclassid,0)
    plhs(2) = mxCreateNumericArray(3,[e5,e6,e7],myclassid,0)
    plhs(3) = mxCreateNumericArray(3,[e5,e6,e7],myclassid,0)
    plhs(4) = mxCreateNumericArray(3,[e5,e6,e7],myclassid,0)
    o1_pr = mxGetPr(plhs(1))      
    o2_pr = mxGetPr(plhs(2))
    o3_pr = mxGetPr(plhs(3))
    o4_pr = mxGetPr(plhs(4))
                      
    ! Call subroutine for assignment
    call allterp(	o1,o2,o3,o4, &
					nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
					x1,x2,x3,x4,x5,x6,x7, &
					x1i,x2i,x3i,x4i,x5i,x6i,x7i, &
					pf1,pf2,pf3,pf4, &
					e5,e6,e7)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(o1,o1_pr,e5*e6*e7)
    call mxCopyReal8ToPtr(o2,o2_pr,e5*e6*e7)
    call mxCopyReal8ToPtr(o3,o3_pr,e5*e6*e7)
    call mxCopyReal8ToPtr(o4,o4_pr,e5*e6*e7)
    
    ! Deallocate arrays
    deallocate(x1)
    deallocate(x2)  
    deallocate(x3) 
    deallocate(x4) 
    deallocate(x5) 
    deallocate(x6) 
    deallocate(x7)
    deallocate(x5i)     
    deallocate(x6i)
    deallocate(x7i)      
    deallocate(pf1)
    deallocate(pf2)
    deallocate(pf3)
    deallocate(pf4)
    deallocate(o1)   
    deallocate(o2)
    deallocate(o3)
    deallocate(o4)
    
end subroutine mexFunction

subroutine allterp(	o1,o2,o3,o4, &
					nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
					x1,x2,x3,x4,x5,x6,x7, &
					x1i,x2i,x3i,x4i,x5i,x6i,x7i, &
					pf1,pf2,pf3,pf4, &
					e5,e6,e7)

    implicit none
    mwSize :: nx1,nx2,nx3,nx4,nx5,nx6,nx7,e5,e6,e7
    mwSize i5,i6,i7
    double precision x1i,x2i,x3i,x4i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4
    double precision, dimension(nx5) :: x5
    double precision, dimension(nx6) :: x6
    double precision, dimension(nx7) :: x7
    double precision, dimension(nx1,nx2,nx3,nx4,nx5,nx6,nx7) :: pf1,pf2,pf3,pf4
    double precision, dimension(e5) :: x5i
    double precision, dimension(e6) :: x6i
    double precision, dimension(e7) :: x7i
    double precision, dimension(e5,e6,e7) :: o1,o2,o3,o4

    ! Interpolate
    do i7 = 1,e7
      do i6 = 1,e6
        do i5 = 1,e5
          o1(i5,i6,i7) = fastrap(nx1,nx2,nx3,nx4,nx5,nx6,nx7,x1,x2,x3,x4,x5,x6,x7,x1i,x2i,x3i,x4i,x5i(i5),x6i(i6),x7i(i7),pf1)
          o2(i5,i6,i7) = fastrap(nx1,nx2,nx3,nx4,nx5,nx6,nx7,x1,x2,x3,x4,x5,x6,x7,x1i,x2i,x3i,x4i,x5i(i5),x6i(i6),x7i(i7),pf2)
          o3(i5,i6,i7) = fastrap(nx1,nx2,nx3,nx4,nx5,nx6,nx7,x1,x2,x3,x4,x5,x6,x7,x1i,x2i,x3i,x4i,x5i(i5),x6i(i6),x7i(i7),pf3)
          o4(i5,i6,i7) = fastrap(nx1,nx2,nx3,nx4,nx5,nx6,nx7,x1,x2,x3,x4,x5,x6,x7,x1i,x2i,x3i,x4i,x5i(i5),x6i(i6),x7i(i7),pf4)
        end do   
      end do
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
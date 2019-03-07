#include "fintrf.h"
!#include "link_fnl_shared.h"
#include "link_fnl_shared_imsl.h"

module structsolution
    
    integer, parameter :: dp = selected_real_kind(15,307)
    
	! Model dimensions
	mwSize :: nstate,ninit,npf,nshock
	! Options
    mwSize :: Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts
	mwSize :: Oe_pts,Ou_pts,Ov_pts
	! Parameters
    real(dp) :: Ptol,Pvarphi,Ppi,Pzlb,Pphipi,Pphiy
    real(dp) :: Ptheta,Prhor,Peta,Ph,Palpha,Pdelta,Pnu,Pg
	! Steady State
    real(dp) :: Sr,Schi
	! Grids
    real(dp), allocatable, dimension(:) :: Gc_grid,Gk_grid,Gi_grid,Grn_grid,Gg_grid,Gbeta_grid,Gmp_grid
    real(dp), allocatable, dimension(:) :: Ge_nodes,Gu_nodes,Gv_nodes
    real(dp), allocatable, dimension(:,:) :: Ge_weight,Gu_weight,Gv_weight
    real(dp), allocatable, dimension(:,:,:,:,:,:,:) :: Gc_gr,Gk_gr,Gi_gr,Grn_gr,Gg_gr,Gmp_gr
	! Policy Functions
	real(dp), allocatable, dimension(:,:,:,:,:,:,:) :: pfn,pfq,pfpi,pfpsi
	
end module structsolution

module csolvepass
    
    use structsolution, only : dp,Oe_pts,Ou_pts,Ov_pts,ninit
    
    real(dp), allocatable, dimension(:,:) :: state
	real(dp), allocatable, dimension(:,:,:) :: gpArr,betapArr
	real(dp), allocatable, dimension(:,:,:) :: weightArr
	
	!$omp threadprivate(state,weightArr)

contains
	subroutine allocateprivate
	
		implicit none
		
		allocate(state(ninit,1))
		allocate(weightArr(Oe_pts,Ou_pts,Ov_pts))
	end subroutine allocateprivate
	
	subroutine deallocateprivate
	
		implicit none
		
		deallocate(state)
		deallocate(weightArr)
	end subroutine deallocateprivate
	
end module csolvepass

subroutine mexFunction(nlhs,plhs,nrhs,prhs)
    ! Declarations
	use structsolution
	implicit none

    ! mexFunction argument
    mwPointer plhs(*),prhs(*)    
    integer*4 nlhs,nrhs      
    
    ! Function declarations
    mwSize mxGetN,mxGetM
    mwpointer mxGetPr,mxCreateNumericArray
    real(dp) mxGetScalar  
    integer*4 mxClassIDFromClassName  
    
    ! Pointers to input
    mwpointer Gc_grid_pr,Gk_grid_pr,Gi_grid_pr,Grn_grid_pr,Gg_grid_pr,Gbeta_grid_pr,Gmp_grid_pr
    mwpointer Ge_weight_pr,Ge_nodes_pr
    mwpointer Gu_weight_pr,Gu_nodes_pr
    mwpointer Gv_weight_pr,Gv_nodes_pr
    mwpointer pfn_pr,pfq_pr,pfpi_pr,pfpsi_pr
    mwpointer Gc_gr_pr,Gk_gr_pr,Gi_gr_pr,Grn_gr_pr,Gg_gr_pr,Gmp_gr_pr

    ! Internal
    mwSize :: nodes
    integer*4 myclassid,myclassid2,myclassid3
    
    ! Output
	mwpointer :: converged_pr,reason_pr,it_pr,dist_last_pr
	mwSize :: converged,reason,it
	real(dp) :: dist_last

    ! Load Inputs
    ! Parameters
    nstate = mxGetScalar(prhs(1))
    ninit = nstate - 1;
    npf = mxGetScalar(prhs(2))
	nshock = mxGetScalar(prhs(3))
    
    Ptol = mxGetScalar(prhs(4))
    Pvarphi = mxGetScalar(prhs(5))
    Ppi = mxGetScalar(prhs(6))
    Pzlb = mxGetScalar(prhs(7))
    Pphipi = mxGetScalar(prhs(8))
    Pphiy = mxGetScalar(prhs(9))
    Ptheta = mxGetScalar(prhs(10))
    Prhor = mxGetScalar(prhs(11))
    Peta = mxGetScalar(prhs(12))
    Ph = mxGetScalar(prhs(13))
    Palpha = mxGetScalar(prhs(14))
    Pdelta = mxGetScalar(prhs(15))
    Pnu = mxGetScalar(prhs(16))
	Pg = mxGetScalar(prhs(17))
    
	Sr = mxGetScalar(prhs(18))
	Schi = mxGetScalar(prhs(19))
    
    ! Grids
    Oc_pts = mxGetN(prhs(20))
    allocate(Gc_grid(Oc_pts))
    Gc_grid_pr = mxGetPr(prhs(20))
    call mxCopyPtrToReal8(Gc_grid_pr,Gc_grid,Oc_pts)
	
    Ok_pts = mxGetN(prhs(21))
    allocate(Gk_grid(Ok_pts))
    Gk_grid_pr = mxGetPr(prhs(21))
    call mxCopyPtrToReal8(Gk_grid_pr,Gk_grid,Ok_pts)
    
    Oi_pts = mxGetN(prhs(22))
    allocate(Gi_grid(Oi_pts))
    Gi_grid_pr = mxGetPr(prhs(22))
    call mxCopyPtrToReal8(Gi_grid_pr,Gi_grid,Oi_pts)    
    
    Orn_pts = mxGetN(prhs(23))
    allocate(Grn_grid(Orn_pts))
    Grn_grid_pr = mxGetPr(prhs(23))
    call mxCopyPtrToReal8(Grn_grid_pr,Grn_grid,Orn_pts)
    
    Og_pts = mxGetN(prhs(24))
    allocate(Gg_grid(Og_pts))
    Gg_grid_pr = mxGetPr(prhs(24))
    call mxCopyPtrToReal8(Gg_grid_pr,Gg_grid,Og_pts)
    
    Obeta_pts = mxGetN(prhs(25))
    allocate(Gbeta_grid(Obeta_pts))
    Gbeta_grid_pr = mxGetPr(prhs(25))
    call mxCopyPtrToReal8(Gbeta_grid_pr,Gbeta_grid,Obeta_pts)
    
    Omp_pts = mxGetN(prhs(26))
    allocate(Gmp_grid(Omp_pts))
    Gmp_grid_pr = mxGetPr(prhs(26))
    call mxCopyPtrToReal8(Gmp_grid_pr,Gmp_grid,Omp_pts)  
    
	! Shocks
    Oe_pts = mxGetM(prhs(27))
	allocate(Ge_weight(Oe_pts,Oe_pts))
    Ge_weight_pr = mxGetPr(prhs(27))
    call mxCopyPtrToReal8(Ge_weight_pr,Ge_weight,Oe_pts*Oe_pts)
    
    allocate(Ge_nodes(Oe_pts))
    Ge_nodes_pr = mxGetPr(prhs(28))
    call mxCopyPtrToReal8(Ge_nodes_pr,Ge_nodes,Oe_pts)
    
    Ou_pts = mxGetM(prhs(29))
	allocate(Gu_weight(Ou_pts,Ou_pts))
    Gu_weight_pr = mxGetPr(prhs(29))
    call mxCopyPtrToReal8(Gu_weight_pr,Gu_weight,Ou_pts*Ou_pts)
    
    allocate(Gu_nodes(Ou_pts))
    Gu_nodes_pr = mxGetPr(prhs(30))
    call mxCopyPtrToReal8(Gu_nodes_pr,Gu_nodes,Ou_pts)    
    
    Ov_pts = mxGetM(prhs(31))
	allocate(Gv_weight(Ov_pts,Ov_pts))
    Gv_weight_pr = mxGetPr(prhs(31))
    call mxCopyPtrToReal8(Gv_weight_pr,Gv_weight,Ov_pts*Ov_pts)
    
    allocate(Gv_nodes(Ov_pts))
    Gv_nodes_pr = mxGetPr(prhs(32))
    call mxCopyPtrToReal8(Gv_nodes_pr,Gv_nodes,Ov_pts)     
    
    ! Policy functions
    nodes = Oc_pts*Ok_pts*Oi_pts*Orn_pts*Og_pts*Obeta_pts*Omp_pts
    allocate(pfn(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(pfq(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(pfpi(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(pfpsi(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(Gc_gr(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts)) 
    allocate(Gk_gr(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts)) 
    allocate(Gi_gr(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts)) 
    allocate(Grn_gr(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(Gg_gr(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(Gmp_gr(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    pfn_pr = mxGetPr(prhs(33))
    pfq_pr = mxGetPr(prhs(34))
    pfpi_pr = mxGetPr(prhs(35))
    pfpsi_pr = mxGetPr(prhs(36))
    Gc_gr_pr = mxGetPr(prhs(37))   
    Gk_gr_pr = mxGetPr(prhs(38))    
    Gi_gr_pr = mxGetPr(prhs(39))  
    Grn_gr_pr = mxGetPr(prhs(40))   
    Gg_gr_pr = mxGetPr(prhs(41))   
    Gmp_gr_pr = mxGetPr(prhs(42)) 
    call mxCopyPtrToReal8(pfn_pr,pfn,nodes) 
    call mxCopyPtrToReal8(pfq_pr,pfq,nodes) 
    call mxCopyPtrToReal8(pfpi_pr,pfpi,nodes) 
    call mxCopyPtrToReal8(pfpsi_pr,pfpsi,nodes)
    call mxCopyPtrToReal8(Gc_gr_pr,Gc_gr,nodes)
    call mxCopyPtrToReal8(Gk_gr_pr,Gk_gr,nodes)
    call mxCopyPtrToReal8(Gi_gr_pr,Gi_gr,nodes)
    call mxCopyPtrToReal8(Grn_gr_pr,Grn_gr,nodes)
    call mxCopyPtrToReal8(Gg_gr_pr,Gg_gr,nodes)
    call mxCopyPtrToReal8(Gmp_gr_pr,Gmp_gr,nodes)
    
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    myclassid2 = mxClassIDFromClassName('int8')
    myclassid3 = mxClassIDFromClassName('int32')
    plhs(1) = mxCreateNumericArray(nstate,[Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts],myclassid,0)
    plhs(2) = mxCreateNumericArray(nstate,[Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts],myclassid,0)
    plhs(3) = mxCreateNumericArray(nstate,[Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts],myclassid,0)
    plhs(4) = mxCreateNumericArray(nstate,[Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts],myclassid,0)
	plhs(5) = mxCreateNumericArray(1,1,myclassid2,0)
	plhs(6) = mxCreateNumericArray(1,1,myclassid2,0)
	plhs(7) = mxCreateNumericArray(1,1,myclassid3,0)
	plhs(8) = mxCreateNumericArray(1,1,myclassid,0)
    pfn_pr = mxGetPr(plhs(1))
    pfq_pr = mxGetPr(plhs(2))
    pfpi_pr = mxGetPr(plhs(3))
    pfpsi_pr = mxGetPr(plhs(4))
    converged_pr = mxGetPr(plhs(5))
    reason_pr = mxGetPr(plhs(6))
    it_pr = mxGetPr(plhs(7))
    dist_last_pr = mxGetPr(plhs(8))

    ! Run script
    call script(converged,reason,it,dist_last)    
    
    ! Load Fortran array to pointer (output to MATLAB)   
    call mxCopyReal8toPtr(pfn,pfn_pr,nodes)
    call mxCopyReal8toPtr(pfq,pfq_pr,nodes)
    call mxCopyReal8toPtr(pfpi,pfpi_pr,nodes)
    call mxCopyReal8toPtr(pfpsi,pfpsi_pr,nodes)
    call mxCopyInteger1ToPtr(converged,converged_pr,1)
    call mxCopyInteger1ToPtr(reason,reason_pr,1)
    call mxCopyInteger1ToPtr(it,it_pr,1)
    call mxCopyReal8toPtr(dist_last,dist_last_pr,1)
    
    ! Deallocate arrays
    deallocate(Gc_grid)
    deallocate(Gk_grid)
    deallocate(Gi_grid)
    deallocate(Grn_grid)
    deallocate(Gg_grid)
    deallocate(Gbeta_grid)
    deallocate(Gmp_grid)
    deallocate(Ge_weight)  
    deallocate(Ge_nodes) 
    deallocate(Gu_weight)  
    deallocate(Gu_nodes) 
    deallocate(Gv_weight)  
    deallocate(Gv_nodes)   
    deallocate(pfn)
    deallocate(pfq)
    deallocate(pfpi)
    deallocate(pfpsi)
    deallocate(Gc_gr)
    deallocate(Gk_gr)
    deallocate(Gi_gr)
    deallocate(Grn_gr)
    deallocate(Gg_gr)
    deallocate(Gmp_gr)
    
end subroutine mexFunction

subroutine script(converged,reason,it,dist_last)
        
    ! Declarations
	use structsolution
	use csolvepass
    implicit none

    ! Outputs
    mwSize, intent(out) :: converged,reason,it
	real(dp), intent(out) :: dist_last
    
    ! Internal
    mwSize i1,i2,i3,i4,i5,i6,i7,j1,j2,j3
    real(dp), allocatable, dimension(:,:,:,:) :: e_weightArr4,u_weightArr4,v_weightArr4
    real(dp), allocatable, dimension(:,:,:,:,:,:,:) :: pfn_up,pfq_up,pfpi_up,pfpsi_up
    real(dp), allocatable, dimension(:,:,:,:,:,:,:) :: dist_n,dist_q,dist_pi,dist_psi
    real(dp) dist_max
    real(dp), dimension(npf,1) :: start,argzero
	real(dp), allocatable, dimension(:,:,:,:,:,:,:) :: y,rgdp,rgdp_gr,rnp,izlb
	real(dp) :: locs,perbind

    ! Allocate memory    
    allocate(e_weightArr4(Oe_pts,Ou_pts,Ov_pts,Og_pts))
    allocate(u_weightArr4(Oe_pts,Ou_pts,Ov_pts,Obeta_pts))
    allocate(v_weightArr4(Oe_pts,Ou_pts,Ov_pts,Omp_pts))
	allocate(pfn_up(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
	allocate(pfq_up(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
	allocate(pfpi_up(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(pfpsi_up(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(dist_n(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))  
	allocate(dist_q(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))  
    allocate(dist_pi(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(dist_psi(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(y(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(rgdp(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
    allocate(rgdp_gr(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
	allocate(rnp(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
	allocate(izlb(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts))
	allocate(gpArr(Oe_pts,Ou_pts,Ov_pts))
	allocate(betapArr(Oe_pts,Ou_pts,Ov_pts))
    
    ! Calculate processes and weights
    do j3 = 1,Ov_pts  
	  do j2 = 1,Ou_pts
        do j1 = 1,Oe_pts
          gpArr(j1,j2,j3) = Ge_nodes(j1);
          betapArr(j1,j2,j3) = Gu_nodes(j2);
        end do
      end do 
    end do
    
	do i5 = 1,Og_pts
      do j3 = 1,Ov_pts  
	    do j2 = 1,Ou_pts
          do j1 = 1,Oe_pts
            e_weightArr4(j1,j2,j3,i5) = Ge_weight(i5,j1);
          end do
        end do 
      end do
    end do
	
	do i6 = 1,Obeta_pts
      do j3 = 1,Ov_pts  
	    do j2 = 1,Ou_pts
          do j1 = 1,Oe_pts
            u_weightArr4(j1,j2,j3,i6) = Gu_weight(i6,j2);
          end do
        end do 
      end do
    end do
	
	do i7 = 1,Omp_pts
      do j3 = 1,Ov_pts  
	    do j2 = 1,Ou_pts
          do j1 = 1,Oe_pts
            v_weightArr4(j1,j2,j3,i7) = Gv_weight(i7,j3);
          end do
        end do 
      end do
    end do   
    
	! Policy Function Iteration Algorithm  
    it = 1;                                 ! Iteration Counter
    converged = -1;                         ! Convergence Flag
	reason = 0; 							! Stopping reason
    ! Time iteration/linear interpolation algorithm
    do while (converged == -1)
        !$omp parallel private(start,argzero)
	    call allocateprivate
        !$omp do collapse(7)
        do i7 = 1,Omp_pts
		  do i6 = 1,Obeta_pts
            do i5 = 1,Og_pts  
              do i4 = 1,Orn_pts
                do i3 = 1,Oi_pts
                  do i2 = 1,Ok_pts
			        do i1 = 1,Oc_pts
			          ! State variables
			          start(1,1) = pfn(i1,i2,i3,i4,i5,i6,i7)
			          start(2,1) = pfq(i1,i2,i3,i4,i5,i6,i7)
                      start(3,1) = pfpi(i1,i2,i3,i4,i5,i6,i7)
                      start(4,1) = pfpsi(i1,i2,i3,i4,i5,i6,i7)
                      state(1,1) = Gc_grid(i1)
				      state(2,1) = Gk_grid(i2)
                      state(3,1) = Gi_grid(i3)
                      state(4,1) = Grn_grid(i4)
                      state(5,1) = Gg_grid(i5) 
                      state(6,1) = Gmp_grid(i7)     
                      weightArr = e_weightArr4(:,:,:,i5)*u_weightArr4(:,:,:,i6)*v_weightArr4(:,:,:,i7)
                      ! Find optimal policy functions on each node.
			          ! csolve finds the zeros of 'eqm'
			          ! Start csolve with the current policy function 
			          argzero = start
			          call csolve(argzero)
				
			          ! Store updated policy functions            
			          pfn_up(i1,i2,i3,i4,i5,i6,i7) = argzero(1,1)
			          pfq_up(i1,i2,i3,i4,i5,i6,i7) = argzero(2,1)
                      pfpi_up(i1,i2,i3,i4,i5,i6,i7) = argzero(3,1)
                      pfpsi_up(i1,i2,i3,i4,i5,i6,i7) = argzero(4,1)
                    end do
                  end do
                end do                
              end do
            end do  
          end do
        end do
        !$omp end do
		call deallocateprivate
        !$omp end parallel
  
        ! Policy function distances
        dist_n = abs(pfn_up - pfn);
        dist_q = abs(pfq_up - pfq);
        dist_pi = abs(pfpi_up - pfpi);
        dist_psi = abs(pfpsi_up - pfpsi);
        
        ! Maximum distance
		dist_max = max(maxval(dist_n),maxval(dist_q),maxval(dist_pi),maxval(dist_psi));
			
        ! Update policy functions
        pfn = pfn_up;
        pfq = pfq_up;
		pfpi = pfpi_up;
        pfpsi = pfpsi_up;
        
		! Stopping reasons
		if (dist_max > 0.1) then
			reason = 1;
		else if (all(pfpi < 0.5) .OR. any(pfn < 0)) then
			reason = 2;
		end if
		
        ! Check convergence criterion
        if (dist_max < Ptol) then 
            converged = 1;
		else if (reason > 0) then
			converged = 0;
		end if
        
        ! Iteration Information
        if (mod(it,10) == 1 .OR. converged >= 0) then
            ! Find where ZLB binds 
            !   Production function
            y = (Gk_gr/Gg_gr)**Palpha*pfn**(1-Palpha);
            ! Real GDP
            rgdp = (1-(Pvarphi*(pfpi/Ppi-1)**2)/2)*y;
            rgdp_gr = Gc_gr + Gi_gr;
            !   Notional interest rate
            rnp = Grn_gr**Prhor*(Sr*(pfpi/Ppi)**Pphipi*(Gg_gr*rgdp/(Pg*rgdp_gr))**Pphiy)**(1-Prhor)*exp(Gmp_gr)
            !   Nodes where ZLB binds
		    izlb = 0.d0
		    where (rnp <= Pzlb) izlb = 1.d0;
            !   Percent nodes binding
            perbind = 100*sum(izlb)/(Oc_pts*Ok_pts*Oi_pts*Orn_pts*Og_pts*Obeta_pts*Omp_pts);
            ! Display iteration information
            call itinfo(it,dist_max,perbind);
        end if
		it = it + 1;
    end do

	it = it - 1;
	dist_last = dist_max;
   
    ! Deallocate memory 
    deallocate(e_weightArr4)
    deallocate(u_weightArr4)
    deallocate(v_weightArr4)
	deallocate(pfn_up)
	deallocate(pfq_up)
    deallocate(pfpi_up)
    deallocate(pfpsi_up)
	deallocate(dist_n)  
	deallocate(dist_q)  
    deallocate(dist_pi)
    deallocate(dist_psi)
    deallocate(y)
    deallocate(rgdp)
    deallocate(rgdp_gr)
    deallocate(rnp)
    deallocate(izlb)
	deallocate(gpArr)
    deallocate(betapArr)
	
end subroutine script

! Root finder based on Chris Sims csolve.m, no verbose
subroutine csolve(x)
    
    use structsolution
	use rnnoa_int
	use rnset_int
    use cond_int
    use eye_int
    use lin_sol_gen_int
    implicit none
    
    ! Inputs
    real(dp), dimension(npf,1), intent(inout) :: x
    
    ! Internal
    real(dp), dimension(npf) :: xtemp
    complex*16, dimension(npf) :: fxtemp
    complex*16, dimension(npf,1) :: fx,fmin
    real(dp), dimension(npf,1) :: fx2,xmin,randn1,dx0,dx
	real(dp), dimension(npf) :: randn0
    real(dp), parameter :: delta = 1e-6, alpha = 1e-3, crit = 1e-4
    mwSize, parameter :: itmax = 10
    real(dp), dimension(npf,npf) :: tvec,guessh,gradinv,zeros
    complex*16, dimension(npf,npf) :: fguessh,grad
    real(dp) :: af0,af00,afmin,af
    mwSize j,itct,shrink,done,subDone
    real(dp) :: lambda,lambdamin,factor,dxSize
	mwSize :: iseed,randomize
    real(dp), dimension(npf,0) :: b0,x0
	real(dp) :: det0(2)
    
    ! Output
    real(dp), dimension(npf,1) :: argzero
	
	iseed = 123456
	call rnset(iseed)
	zeros = 0.0d0;
	randomize = 0;
    tvec = delta*eye(npf);
    done = 0;
    call eqm(x,1,fx);
    af0 = sum(abs(fx));
    af00 = af0;
    itct = 0;
    do while (done == 0)
        xtemp(1) = x(1,1)
        xtemp(2) = x(2,1)
        xtemp(3) = x(3,1)
        xtemp(4) = x(4,1)
        guessh = spread(xtemp,2,npf)+tvec
		if (itct>3 .AND. af00-af0<crit*max(1.d0,af0) .AND. mod(itct,2)==1) then
			randomize = 1;
		else
			call eqm(guessh,npf,fguessh);
			fxtemp(1) = fx(1,1)
			fxtemp(2) = fx(2,1)
            fxtemp(3) = fx(3,1)
            fxtemp(4) = fx(4,1)
			grad = (fguessh - spread(fxtemp,2,npf))/delta;
			if (.FALSE. .EQV. any(aimag(grad) .NE. zeros)) then
				call D_LIN_SOL_GEN(real(grad),b0,x0,nrhs=0,ainv=gradinv,det=det0)
				if (abs(det0(1))**det0(2) < 1.0D-12) then
					grad = grad+tvec;
                    call D_LIN_SOL_GEN(real(grad),b0,x0,nrhs=0,ainv=gradinv,det=det0)
				end if
				dx0 = matmul(-gradinv,fx);
                randomize=0;
			else
				randomize=1;
			end if
		end if
	    if (randomize == 1) then
			call rnnoa(randn0)
			randn1(1,1) = randn0(1)
			randn1(2,1) = randn0(2)
            randn1(3,1) = randn0(3)
            randn1(4,1) = randn0(4)
			dx0=sqrt(sum(x**2))/randn1;
	    end if
        lambda = 1;
        lambdamin = 1;
        fmin = fx;
        xmin = x;
        afmin = af0;
        dxSize = sqrt(sum(dx0**2))
        factor = 0.6d0;
        shrink = 1;
        subDone = 0;
        do while (subDone == 0)
            dx = lambda*dx0;
            call eqm(x+dx,1,fx);
            af = sum(abs(fx));
            if (af < afmin) then
                afmin = af;
                fmin = fx;
                xmin = x + dx;
            end if
          if (((lambda > 0) .AND. (af0-af < alpha*lambda*af0)) .OR. ((lambda < 0) .AND. (af0-af < 0))) then
             if (shrink == 0) then
                factor = factor**0.6d0;
                shrink = 1;
             end if
             if (abs(lambda*(1-factor))*dxSize > .1*delta) then
                lambda = factor*lambda;
             else if ((lambda > 0) .AND. (factor==0.6d0)) then
                lambda = -0.3d0;
             else
                subDone = 1;
             end if
          else if ((lambda >0) .AND. (af-af0 > (1-alpha)*lambda*af0)) then
             if (shrink == 1) then
                factor = factor**0.6d0;
                shrink = 0;
             end if
             lambda = lambda/factor;
          else
             subDone=1;
          end if
       end do
       itct = itct+1;
       x = xmin;
       fx = fmin;
       af00 = af0;
       af0 = afmin;
       if ((itct >= itmax) .OR. (af0 < crit)) then
          done = 1;
       end if
    end do
    
    argzero = x

end subroutine csolve

subroutine eqm(x,ncol,Res)
   
	use structsolution
	use csolvepass
    implicit none
	
    ! Inputs
    mwSize, intent(in) :: ncol
    real(dp), dimension(npf,ncol), intent(in) :: x
    
    ! Output
    complex*16, dimension(npf,ncol) :: Res
    
    ! Internal
    real(dp) :: c,k,i,rn,g,mp,n,q,pie,psi,rgdp,r
    complex*16 :: y,rgdpp,rnp,w,cp,ip,iac,kp
	real(dp), dimension(Oe_pts,Ou_pts,Ov_pts) :: npArr,qpArr,pipArr,psipArr
    complex*16, dimension(Oe_pts,Ou_pts,Ov_pts) :: ypArr,rkpArr,wpArr,cppArr,ippArr,iacpArr,sdfArr
    real(dp) :: Ecap,Einv,Ebond,Efp
    mwSize :: j

	! State Values	
    c = state(1,1);   	  !Consumption last period
    k = state(2,1);  	  !Capital last period
    i = state(3,1);       !Investment last period   
    rn = state(4,1);  	  !Notional interest rate last period
    g = state(5,1);       !Growth state current period
    mp  = state(6,1);     !Monetary policy shock current state
    
	do j = 1,ncol  
		! Policy Function Guesses
        n = x(1,j);     !Labor policy current period 
        q = x(2,j);     !Tobin's q current period
        pie = x(3,j);   !Inflation current period  
        psi = x(4,j);   !Marginal cost current period
		!----------------------------------------------------------------------
		! Solve for variables
		!----------------------------------------------------------------------
		! Production function
        y = (k/g)**Palpha*n**(1-Palpha); 
        ! Real GDP        
        rgdp = c+i;
        rgdpp = (1-Pvarphi*(pie/Ppi-1)**2/2)*y;
        ! Interest rate rule
        rnp = rn**Prhor*(Sr*(pie/Ppi)**Pphipi*(g*rgdpp/(Pg*rgdp))**Pphiy)**(1-Prhor)*exp(mp);
		r = max(Pzlb,abs(rnp));
        ! Firm FOC labor
        w = (1-Palpha)*psi*y/n;
        ! FOC labor
        cp = w/(Schi*n**Peta)+Ph*c/g;
        ! Aggregate resource constraint  
        ip = rgdpp - cp; 
        ! Investment adjustment costs
        iac = ip*g/(i*Pg);
        ! Law of motion for capital
        kp = (1-Pdelta)*(k/g)+ip*(1-Pnu*(iac-1)**2/2);
		!----------------------------------------------------------------------
		! Linear interpolation of the policy variables 
		!----------------------------------------------------------------------
        call allterp743(Oe_pts,Ou_pts,Ov_pts, &
                        Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts,Omp_pts, &
                        Gc_grid,Gk_grid,Gi_grid,Grn_grid, &
                        cp,kp,ip,rnp, &
                        pfn,pfq,pfpi,pfpsi, & 
                        npArr,qpArr,pipArr,psipArr)           
		!----------------------------------------------------------------------      
        ! Solve for variables inside expectations
		!----------------------------------------------------------------------
	    ! Production function
        ypArr = (kp/gpArr)**Palpha*npArr**(1-Palpha);
        ! Firm FOC capital
        rkpArr = Palpha*psipArr*gpArr*ypArr/kp;
        ! Firm FOC labor
        wpArr = (1-Palpha)*psipArr*ypArr/npArr;
        ! FOC labor
        cppArr = wpArr/(Schi*npArr**Peta)+Ph*cp/gpArr;
        ! Aggregate resource constraint  
        ippArr = ypArr*(1-Pvarphi*(pipArr/Ppi-1)**2/2) - cppArr;   
        ! Investment adjustment costs
        iacpArr = ippArr*gpArr/(ip*Pg);
        ! Stochastic discount factor
        sdfArr = betapArr*(cp-Ph*c/g)/(cppArr-Ph*cp/gpArr);
		!----------------------------------------------------------------------
		! Numerical integration
		!----------------------------------------------------------------------
        ! Compute all combinations of shocks
        Ecap = sum(weightArr*(sdfArr/gpArr)*(rkpArr+(1-Pdelta)*qpArr));
        Einv = sum(weightArr*sdfArr*qpArr*(gpArr/Pg)*(ippArr/ip)**2*(iacpArr-1));   
        Ebond = sum(weightArr*sdfArr/(gpArr*pipArr));
        Efp = sum(weightArr*sdfArr*(pipArr/Ppi-1)*(ypArr/y)*(pipArr/Ppi));     
		!----------------------------------------------------------------------
		! First-order conditions
		!----------------------------------------------------------------------
		! FOC capital
        Res(1,j) = q - Ecap;
        ! FOC investment
        Res(2,j) = 1-q*(1-Pnu*(iac-1)**2/2-Pnu*iac*(iac-1))-Pnu*Einv;    
        ! FOC bond
        Res(3,j) = 1 - r*Ebond;
        ! Firm Pricing
        Res(4,j) = Pvarphi*(pie/Ppi-1)*pie/Ppi-(1-Ptheta)-Ptheta*psi-Pvarphi*Efp;
    end do

end subroutine eqm  
 
subroutine allterp743(Oe_pts,Ou_pts,Ov_pts, &
                      nx1,nx2,nx3,nx4,nx5,nx6,nx7, &
                      x1,x2,x3,x4, &
                      x1i,x2i,x3i,x4i, &
                      z1,z2,z3,z4, &
                      o1,o2,o3,o4)
    
    use structsolution, only : dp
    implicit none
	
	! Inputs
    mwSize, intent(in) :: Oe_pts,Ou_pts,Ov_pts,nx1,nx2,nx3,nx4,nx5,nx6,nx7
    real(dp), intent(in) :: x1i,x2i,x3i,x4i
    real(dp), intent(in) :: x1(nx1),x2(nx2),x3(nx3),x4(nx4)
    real(dp), dimension(nx1,nx2,nx3,nx4,nx5,nx6,nx7), intent(in) :: z1,z2,z3,z4
	
	! Outputs
    real(dp), dimension(Oe_pts,Ou_pts,Ov_pts), intent(out) :: o1,o2,o3,o4
    
	! Internal
    real(dp) :: s1, s2, s3, s4
    real(dp) :: x1i_min, x2i_min, x3i_min, x4i_min
    mwSize loc1, loc2, loc3, loc4
    real(dp), dimension(4) :: xi, xi_left, xi_right, w_2, w_1
    real(dp), dimension(2) :: w1, w2, w3, w4
    mwSize m1, m2, m3, m4, i5, i6, i7
	real(dp) :: wtemp
    
    ! Grid intervals
    s1 = x1(2) - x1(1)
    s2 = x2(2) - x2(1)
    s3 = x3(2) - x3(1)
    s4 = x4(2) - x4(1)
    
    ! Endogenous state variables   
    x1i_min = x1i - x1(1)
    loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));
    
    x2i_min = x2i - x2(1)
    loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));
    
    x3i_min = x3i - x3(1)
    loc3 = min(nx3-1,max(1,floor(x3i_min/s3) + 1));
    
    x4i_min = x4i - x4(1)
    loc4 = min(nx4-1,max(1,floor(x4i_min/s4) + 1));
    
    xi = [x1i, x2i, x3i, x4i]
    xi_left = [x1(loc1), x2(loc2), x3(loc3), x4(loc4)]
    xi_right = [x1(loc1+1), x2(loc2+1), x3(loc3+1), x4(loc4+1)]
    
    w_2 = (xi - xi_left)/(xi_right - xi_left)
    w_1 = 1 - w_2
    w1 = [w_1(1), w_2(1)]
    w2 = [w_1(2), w_2(2)]
    w3 = [w_1(3), w_2(3)]
    w4 = [w_1(4), w_2(4)]
    
    o1 = 0.d0
    o2 = 0.d0
    o3 = 0.d0
    o4 = 0.d0
    do i7 = 1,Ov_pts
      do i6 = 1,Ou_pts
        do i5 = 1,Oe_pts 
          do m4 = 0,1 
            do m3 = 0,1 
              do m2 = 0,1    
                do m1 = 0,1
                  wtemp = w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)  
                  o1(i5,i6,i7) = o1(i5,i6,i7) + wtemp*z1(loc1+m1,loc2+m2,loc3+m3,loc4+m4,i5,i6,i7)
                  o2(i5,i6,i7) = o2(i5,i6,i7) + wtemp*z2(loc1+m1,loc2+m2,loc3+m3,loc4+m4,i5,i6,i7)
                  o3(i5,i6,i7) = o3(i5,i6,i7) + wtemp*z3(loc1+m1,loc2+m2,loc3+m3,loc4+m4,i5,i6,i7)
                  o4(i5,i6,i7) = o4(i5,i6,i7) + wtemp*z4(loc1+m1,loc2+m2,loc3+m3,loc4+m4,i5,i6,i7)
                end do  
              end do
            end do            
          end do
        end do
      end do
    end do
	 
end subroutine allterp743     
    
! Calculates algorithm and iteration duration    
subroutine itinfo(it,dist_max,perbind)

    use structsolution, only : dp
    implicit none    
    ! Inputs
    mwSize it
    real(dp) :: dist_max,perbind
    
    ! Text Display
    integer*4, external :: mexprintf,mexEvalString
    integer*4 pm,evalout
    character*100 line
    
    ! Output iteration information
    write(line,'(A,I6,A,E11.4,A,F7.4,A)') &
        'iter: ', it, &
        ' dist: ', dist_max, &
        ' bind: ', perbind, '%'
    pm = mexPrintf(line//achar(13)) 
    evalout = mexEvalString("pause(.0001);")

end subroutine itinfo
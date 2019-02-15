#include "fintrf.h"
#include "link_fnl_shared.h"
!#include "link_fnl_shared_imsl.h"

module structsimulation
    
    integer, parameter :: dp = selected_real_kind(15,307)
    
    ! Simulation parameters, initial vector, and shocks
    mwSize :: npers,Vnplotvar,nsims
	real(dp), allocatable, dimension(:) :: initvec
	real(dp), allocatable, dimension(:,:) :: e,u
	! Options
    mwSize :: Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts
	! Parameters
    real(dp) :: Pvarphi,Ppi,Pzlb,Pphipi,Pphiy,Ptheta,Prhor
    real(dp) :: Peta,Ph,Palpha,Pdelta,Pnu,Pg,Prhog,Pbeta,Prhobeta   
	! Steady State
    real(dp) :: Sc,Sk,Si,Sr,Schi
    ! Variables
	mwSize :: Vg,Vbeta,Vn,Vpi,Vpsi,Vy,Vrn,Vr,Vw,Vc,Vi,Vk,Vrgdp
	! Grids
    real(dp), allocatable, dimension(:) :: Gc_grid,Gk_grid,Gi_grid,Grn_grid,Gg_grid,Gbeta_grid
	! Policy Functions
	real(dp), allocatable, dimension(:,:,:,:,:,:) :: pfn,pfpi,pfpsi
	
end module structsimulation

subroutine mexFunction(nlhs,plhs,nrhs,prhs)
    ! Declarations
    use structsimulation
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
    mwpointer initvec_pr,Gc_grid_pr,Gk_grid_pr,Gi_grid_pr,Grn_grid_pr,Gg_grid_pr,Gbeta_grid_pr
    mwpointer pfn_pr,pfpi_pr,pfpsi_pr
    mwpointer e_pr,u_pr
					
    ! Internal
    mwSize :: nodes
    integer*4 myclassid
    
    ! Output
    mwpointer :: sims_pr
	real(dp), allocatable, dimension(:,:,:) :: sims

    ! Load Inputs
	
	! Initial state vector
    Vnplotvar = mxGetM(prhs(1))
    allocate(initvec(Vnplotvar))
    initvec_pr = mxGetPr(prhs(1))
    call mxCopyPtrToReal8(initvec_pr,initvec,Vnplotvar)
	
    ! 	Parameters
    Pvarphi = mxGetScalar(prhs(2))
    Ppi = mxGetScalar(prhs(3))
    Pzlb = mxGetScalar(prhs(4))
    Pphipi = mxGetScalar(prhs(5))
    Pphiy = mxGetScalar(prhs(6))
    Ptheta = mxGetScalar(prhs(7))
    Prhor = mxGetScalar(prhs(8))
    Peta = mxGetScalar(prhs(9))
    Ph = mxGetScalar(prhs(10))
    Palpha = mxGetScalar(prhs(11))
    Pdelta = mxGetScalar(prhs(12))
    Pnu = mxGetScalar(prhs(13))
    Pg = mxGetScalar(prhs(14))
    Prhog = mxGetScalar(prhs(15))
    Pbeta = mxGetScalar(prhs(16))
    Prhobeta = mxGetScalar(prhs(17))  
    
    !	Steady State
    Sc = mxGetScalar(prhs(18))
    Sk = mxGetScalar(prhs(19))
    Si = mxGetScalar(prhs(20))
    Sr = mxGetScalar(prhs(21))
    Schi = mxGetScalar(prhs(22))
	
	!	Variables
    Vg = mxGetScalar(prhs(23))
    Vbeta = mxGetScalar(prhs(24))
    Vn = mxGetScalar(prhs(25))
    Vpi = mxGetScalar(prhs(26))
    Vpsi = mxGetScalar(prhs(27))
    Vy = mxGetScalar(prhs(28))
    Vrn = mxGetScalar(prhs(29))
    Vr = mxGetScalar(prhs(30))
    Vw = mxGetScalar(prhs(31))
    Vc = mxGetScalar(prhs(32))
    Vi = mxGetScalar(prhs(33))
    Vk = mxGetScalar(prhs(34))
    Vrgdp = mxGetScalar(prhs(35))
    
    ! 	Grids
    Oc_pts = mxGetN(prhs(36))
    allocate(Gc_grid(Oc_pts))
    Gc_grid_pr = mxGetPr(prhs(36))
    call mxCopyPtrToReal8(Gc_grid_pr,Gc_grid,Oc_pts)
    
    Ok_pts = mxGetN(prhs(37))
    allocate(Gk_grid(Ok_pts))
    Gk_grid_pr = mxGetPr(prhs(37))
    call mxCopyPtrToReal8(Gk_grid_pr,Gk_grid,Ok_pts)
    
    Oi_pts = mxGetN(prhs(38))
    allocate(Gi_grid(Oi_pts))
    Gi_grid_pr = mxGetPr(prhs(38))
    call mxCopyPtrToReal8(Gi_grid_pr,Gi_grid,Oi_pts)        
    
    Orn_pts = mxGetN(prhs(39))
    allocate(Grn_grid(Orn_pts))
    Grn_grid_pr = mxGetPr(prhs(39))
    call mxCopyPtrToReal8(Grn_grid_pr,Grn_grid,Orn_pts)
    
    Og_pts = mxGetN(prhs(40))
    allocate(Gg_grid(Og_pts))
    Gg_grid_pr = mxGetPr(prhs(40))
    call mxCopyPtrToReal8(Gg_grid_pr,Gg_grid,Og_pts)
    
    Obeta_pts = mxGetN(prhs(41))
    allocate(Gbeta_grid(Obeta_pts))
    Gbeta_grid_pr = mxGetPr(prhs(41))
    call mxCopyPtrToReal8(Gbeta_grid_pr,Gbeta_grid,Obeta_pts)
    
    ! Policy functions
    nodes = Oc_pts*Ok_pts*Oi_pts*Orn_pts*Og_pts*Obeta_pts
    allocate(pfn(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts))
    allocate(pfpi(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts))
    allocate(pfpsi(Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts))
	pfn_pr = mxGetPr(prhs(42))
    pfpi_pr = mxGetPr(prhs(43))
    pfpsi_pr = mxGetPr(prhs(44))
    call mxCopyPtrToReal8(pfn_pr,pfn,nodes) 
    call mxCopyPtrToReal8(pfpi_pr,pfpi,nodes) 
    call mxCopyPtrToReal8(pfpsi_pr,pfpsi,nodes)
	
	! 	Shocks
	npers = mxGetM(prhs(45))
	nsims = mxGetN(prhs(45))
    allocate(e(npers,nsims))
    allocate(u(npers,nsims))
	
    e_pr = mxGetPr(prhs(45))
    u_pr = mxGetPr(prhs(46))
    call mxCopyPtrToReal8(e_pr,e,npers*nsims) 
    call mxCopyPtrToReal8(u_pr,u,npers*nsims)
    
    ! Create outputs
    myclassid = mxClassIDFromClassName('double')
    plhs(1) = mxCreateNumericArray(3,[npers,Vnplotvar,nsims],myclassid,0)
    sims_pr = mxGetPr(plhs(1))
    allocate(sims(npers,Vnplotvar,nsims))

    ! Run script
    call simulation(sims)     
    
    ! Load Fortran array to pointer (output to MATLAB)   
    call mxCopyReal8toPtr(sims,sims_pr,npers*Vnplotvar*nsims)
    
    ! Deallocate arrays
    deallocate(initvec)
    deallocate(Gc_grid)
    deallocate(Gk_grid)
    deallocate(Gi_grid)
    deallocate(Grn_grid)
    deallocate(Gg_grid)
    deallocate(Gbeta_grid)
    deallocate(pfn)
    deallocate(pfpi)
    deallocate(pfpsi)
    deallocate(e)
    deallocate(u)
    deallocate(sims)
    
end subroutine mexFunction

subroutine simulation(sims)   
        
    ! Declarations
    use structsimulation
    implicit none

    ! Inputs
	real(dp), dimension(npers,Vnplotvar,nsims), intent(out) :: sims
    
    ! Internal
    mwSize isim,t

	! Initialize state variables at deterministic steady state
	if (all(initvec == 0)) then
		sims(1,Vc,:) = Sc;
        sims(1,Vk,:) = Sk;
        sims(1,Vi,:) = Si;
        sims(1,Vrn,:) = Sr;
        sims(1,Vg,:) = Pg;
		sims(1,Vbeta,:) = Pbeta;
        sims(1,Vrgdp,:) = Sc+Si;
	else
		sims(1,:,:) = spread(initvec,2,nsims);
    end if

	!----------------------------------------------------------------------
	! Simulate time series
	!----------------------------------------------------------------------
	!$omp parallel default(shared)
	!$omp do
	do isim = 1,nsims
		do t = 2,npers
			! Growth rate
			sims(t,Vg,isim) = Pg*(sims(t-1,Vg,isim)/Pg)**Prhog*exp(e(t,isim));
			! Discount factor
			sims(t,Vbeta,isim) = Pbeta*(sims(t-1,Vbeta,isim)/Pbeta)**Prhobeta*exp(u(t,isim));
			! Evaluate policy functions
			call interpolate63(	Oc_pts,Ok_pts,Oi_pts,Orn_pts,Og_pts,Obeta_pts, &
								Gc_grid,Gk_grid,Gi_grid,Grn_grid,Gg_grid,Gbeta_grid, &
								sims(t-1,Vc,isim),sims(t-1,Vk,isim),sims(t-1,Vi,isim), &
                                sims(t-1,Vrn,isim),sims(t,Vg,isim),sims(t,Vbeta,isim), &
								pfn,pfpi,pfpsi, &
								sims(t,Vn,isim),sims(t,Vpi,isim),sims(t,Vpsi,isim));	
			! Production function
            sims(t,Vy,isim) = (sims(t-1,Vk,isim)/sims(t,Vg,isim))**Palpha*sims(t,Vn,isim)**(1-Palpha); 
			! Real GDP
            sims(t,Vrgdp,isim) = (1-Pvarphi*(sims(t,Vpi,isim)/Ppi-1)**2/2)*sims(t,Vy,isim);
            ! Interest rate rule    
			sims(t,Vrn,isim) = sims(t-1,Vrn,isim)**Prhor*(Sr*(sims(t,Vpi,isim)/Ppi)**Pphipi*(sims(t,Vg,isim)*sims(t,Vrgdp,isim)/(Pg*sims(t-1,Vrgdp,isim)))**Pphiy)**(1-Prhor);
			sims(t,Vr,isim) = max(Pzlb,sims(t,Vrn,isim));
            ! Firm FOC labor
            sims(t,Vw,isim) = (1-Palpha)*sims(t,Vpsi,isim)*sims(t,Vy,isim)/sims(t,Vn,isim);
            ! FOC labor
            sims(t,Vc,isim) = sims(t,Vw,isim)/(Schi*sims(t,Vn,isim)**Peta)+Ph*sims(t-1,Vc,isim)/sims(t,Vg,isim);
            ! Aggregate resource constraint
            sims(t,Vi,isim) = sims(t,Vrgdp,isim) - sims(t,Vc,isim);    
            ! Law of motion for capital
            sims(t,Vk,isim) = (1-Pdelta)*(sims(t-1,Vk,isim)/sims(t,Vg,isim))+sims(t,Vi,isim)*(1-Pnu*(sims(t,Vi,isim)*sims(t,Vg,isim)/(Pg*sims(t-1,Vi,isim))-1)**2/2);      
		end do
	end do
	!$omp end do
	!$omp end parallel
	
end subroutine simulation
      
subroutine interpolate63(	nx1,nx2,nx3,nx4,nx5,nx6, &
							x1,x2,x3,x4,x5,x6, &
							x1i,x2i,x3i,x4i,x5i,x6i, &
							z1,z2,z3, &
							o1,o2,o3)

    ! Declarations
    use structsimulation, only : dp
    implicit none
	
	! Inputs
    mwSize, intent(in) :: nx1,nx2,nx3,nx4,nx5,nx6
    real(dp), intent(in) :: x1i,x2i,x3i,x4i,x5i,x6i
    real(dp), intent(in) :: x1(nx1),x2(nx2),x3(nx3),x4(nx4),x5(nx5),x6(nx6)
    real(dp), dimension(nx1,nx2,nx3,nx4,nx5,nx6), intent(in) :: z1,z2,z3
	
	! Outputs
    real(dp), intent(out) :: o1,o2,o3
    
	! Internal
    real(dp) :: s1, s2, s3, s4, s5, s6
    real(dp) :: x1i_min, x2i_min, x3i_min, x4i_min, x5i_min, x6i_min
    mwSize loc1, loc2, loc3, loc4, loc5, loc6
    real(dp), dimension(6) :: xi, xi_left, xi_right, w_2, w_1
    real(dp), dimension(2) :: w1, w2, w3, w4, w5, w6
    mwSize m1, m2, m3, m4, m5, m6, isim
    
    o1 = 0.d0
    o2 = 0.d0
    o3 = 0.d0
	
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

	xi = [x1i, x2i, x3i, x4i, x5i, x6i]
	xi_left = [x1(loc1), x2(loc2), x3(loc3), x4(loc4), x5(loc5), x6(loc6)]
	xi_right = [x1(loc1+1), x2(loc2+1), x3(loc3+1), x4(loc4+1), x5(loc5+1), x6(loc6+1)]

	w_2 = (xi - xi_left)/(xi_right - xi_left)
	w_1 = 1 - w_2
	w1 = [w_1(1), w_2(1)]
	w2 = [w_1(2), w_2(2)]
	w3 = [w_1(3), w_2(3)]
    w4 = [w_1(4), w_2(4)]
	w5 = [w_1(5), w_2(5)]
    w6 = [w_1(6), w_2(6)]
    
    do m6 = 0, 1
      do m5 = 0, 1
        do m4 = 0, 1
          do m3 = 0, 1
            do m2 = 0, 1
              do m1 = 0, 1    
                o1 = o1 + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*w6(m6+1)*z1(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5,loc6+m6) 
	            o2 = o2 + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*w6(m6+1)*z2(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5,loc6+m6) 
	            o3 = o3 + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*w6(m6+1)*z3(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5,loc6+m6)    
              end do
            end do            
          end do
        end do          
      end do
    end do

end subroutine interpolate63
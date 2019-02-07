#include "fintrf.h"
#include "link_fnl_shared_imsl.h"

module structsolution
    
    integer, parameter :: dp = selected_real_kind(15,307)
    
	! Model dimensions
	mwSize :: nstate,npf,nshock
	! Options
    mwSize :: Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts
	mwSize :: Oe_pts,Ou_pts,Ov_pts
	! Parameters
    real(dp) :: Ptol,Pbeta,Pthetap,Peta,Pvarphip,Pphipi,Pphiy,Ph,Prhoi,Prhos,Pg,Ppi,Ps
	! Steady State
    real(dp) :: Schi,Si
	! Grids
    real(dp), allocatable, dimension(:) :: Gg_grid,Gs_grid,Gmp_grid,Gin_grid,Gc_grid
    real(dp), allocatable, dimension(:) :: Ge_nodes,Gu_nodes,Gv_nodes
    real(dp), allocatable, dimension(:,:) :: Ge_weight,Gu_weight,Gv_weight
	real(dp), allocatable, dimension(:,:,:,:,:) :: Gg_gr,Gmp_gr,Gin_gr,Gc_gr
	! Policy Functions
	real(dp), allocatable, dimension(:,:,:,:,:) :: pfc,pfpigap
	
end module structsolution

module csolvepass
    
    use structsolution, only : dp,Oe_pts,Ou_pts,Ov_pts,nstate
    
    real(dp), allocatable, dimension(:,:) :: state
	real(dp), allocatable, dimension(:,:,:) :: gpArr3
	real(dp), allocatable, dimension(:,:,:) :: weightArr3
	
	!$omp threadprivate(state,weightArr3)

contains
	subroutine allocateprivate
	
		implicit none
		
		allocate(state(nstate,1))
		allocate(weightArr3(Oe_pts,Ou_pts,Ov_pts))
	end subroutine allocateprivate
	
	subroutine deallocateprivate
	
		implicit none
		
		deallocate(state)
		deallocate(weightArr3)
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
    mwpointer Gg_grid_pr,Gs_grid_pr,Gmp_grid_pr,Gin_grid_pr,Gc_grid_pr
    mwpointer Ge_weight_pr,Ge_nodes_pr
    mwpointer Gu_weight_pr,Gu_nodes_pr
    mwpointer Gv_weight_pr,Gv_nodes_pr
    mwpointer pfc_pr,pfpigap_pr
    mwpointer Gg_gr_pr,Gmp_gr_pr,Gin_gr_pr,Gc_gr_pr
		
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
    npf = mxGetScalar(prhs(2))
	nshock = mxGetScalar(prhs(3))
    
    Ptol = mxGetScalar(prhs(4))
	Pbeta = mxGetScalar(prhs(5))
	Pthetap = mxGetScalar(prhs(6))
	Peta = mxGetScalar(prhs(7))
	Pg = mxGetScalar(prhs(8))
	Ppi = mxGetScalar(prhs(9))
	Ps = mxGetScalar(prhs(10))
	
	Pvarphip = mxGetScalar(prhs(11))
	Ph = mxGetScalar(prhs(12))
	Prhos = mxGetScalar(prhs(13))
	Prhoi = mxGetScalar(prhs(14))
	Pphipi = mxGetScalar(prhs(15))
	Pphiy = mxGetScalar(prhs(16))
    
	Schi = mxGetScalar(prhs(17))
	Si = mxGetScalar(prhs(18))
	
    ! Grids
    Og_pts = mxGetN(prhs(19))
    allocate(Gg_grid(Og_pts))
    Gg_grid_pr = mxGetPr(prhs(19))
    call mxCopyPtrToReal8(Gg_grid_pr,Gg_grid,Og_pts)
    
    Os_pts = mxGetN(prhs(20))
    allocate(Gs_grid(Os_pts))
    Gs_grid_pr = mxGetPr(prhs(20))
    call mxCopyPtrToReal8(Gs_grid_pr,Gs_grid,Os_pts)
	
    Omp_pts = mxGetN(prhs(21))
    allocate(Gmp_grid(Omp_pts))
    Gmp_grid_pr = mxGetPr(prhs(21))
    call mxCopyPtrToReal8(Gmp_grid_pr,Gmp_grid,Omp_pts) 
	
    Oin_pts = mxGetN(prhs(22))
    allocate(Gin_grid(Oin_pts))
    Gin_grid_pr = mxGetPr(prhs(22))
    call mxCopyPtrToReal8(Gin_grid_pr,Gin_grid,Oin_pts)
	
    Oc_pts = mxGetN(prhs(23))
    allocate(Gc_grid(Oc_pts))
    Gc_grid_pr = mxGetPr(prhs(23))
    call mxCopyPtrToReal8(Gc_grid_pr,Gc_grid,Oc_pts)
    
	! Shocks
    Oe_pts = mxGetM(prhs(24))
	allocate(Ge_weight(Oe_pts,Oe_pts))
    Ge_weight_pr = mxGetPr(prhs(24))
    call mxCopyPtrToReal8(Ge_weight_pr,Ge_weight,Oe_pts*Oe_pts)
    
    allocate(Ge_nodes(Oe_pts))
    Ge_nodes_pr = mxGetPr(prhs(25))
    call mxCopyPtrToReal8(Ge_nodes_pr,Ge_nodes,Oe_pts)
    
    Ou_pts = mxGetM(prhs(26))
	allocate(Gu_weight(Ou_pts,Ou_pts))
    Gu_weight_pr = mxGetPr(prhs(26))
    call mxCopyPtrToReal8(Gu_weight_pr,Gu_weight,Ou_pts*Ou_pts)
    
    allocate(Gu_nodes(Ou_pts))
    Gu_nodes_pr = mxGetPr(prhs(27))
    call mxCopyPtrToReal8(Gu_nodes_pr,Gu_nodes,Ou_pts)
	
    Ov_pts = mxGetM(prhs(28))
	allocate(Gv_weight(Ov_pts,Ov_pts))
    Gv_weight_pr = mxGetPr(prhs(28))
    call mxCopyPtrToReal8(Gv_weight_pr,Gv_weight,Ov_pts*Ov_pts)
    
    allocate(Gv_nodes(Ov_pts))
    Gv_nodes_pr = mxGetPr(prhs(29))
    call mxCopyPtrToReal8(Gv_nodes_pr,Gv_nodes,Ov_pts)      
    
    ! Policy functions
    nodes = Og_pts*Os_pts*Omp_pts*Oin_pts*Oc_pts
    allocate(pfc(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
    allocate(pfpigap(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
    allocate(Gc_gr(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
    allocate(Gin_gr(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))   
    allocate(Gg_gr(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
    allocate(Gmp_gr(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
    pfc_pr = mxGetPr(prhs(30))
    pfpigap_pr = mxGetPr(prhs(31))
    Gg_gr_pr = mxGetPr(prhs(32))  
    Gmp_gr_pr = mxGetPr(prhs(33))
    Gin_gr_pr = mxGetPr(prhs(34))  
    Gc_gr_pr = mxGetPr(prhs(35)) 
    call mxCopyPtrToReal8(pfc_pr,pfc,nodes) 
    call mxCopyPtrToReal8(pfpigap_pr,pfpigap,nodes)
    call mxCopyPtrToReal8(Gg_gr_pr,Gg_gr,nodes)
    call mxCopyPtrToReal8(Gmp_gr_pr,Gmp_gr,nodes)
    call mxCopyPtrToReal8(Gin_gr_pr,Gin_gr,nodes)
    call mxCopyPtrToReal8(Gc_gr_pr,Gc_gr,nodes)
        
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    myclassid2 = mxClassIDFromClassName('int8')
    myclassid3 = mxClassIDFromClassName('int32')
    plhs(1) = mxCreateNumericArray(nstate,[Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts],myclassid,0)
    plhs(2) = mxCreateNumericArray(nstate,[Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts],myclassid,0)
	plhs(3) = mxCreateNumericArray(1,1,myclassid2,0)
	plhs(4) = mxCreateNumericArray(1,1,myclassid2,0)
	plhs(5) = mxCreateNumericArray(1,1,myclassid3,0)
	plhs(6) = mxCreateNumericArray(1,1,myclassid,0)
    pfc_pr = mxGetPr(plhs(1))
    pfpigap_pr = mxGetPr(plhs(2))
    converged_pr = mxGetPr(plhs(3))
    reason_pr = mxGetPr(plhs(4))
    it_pr = mxGetPr(plhs(5))
    dist_last_pr = mxGetPr(plhs(6))

    ! Run script
    call script(converged,reason,it,dist_last)    
    
    ! Load Fortran array to pointer (output to MATLAB)   
    call mxCopyReal8toPtr(pfc,pfc_pr,nodes)
    call mxCopyReal8toPtr(pfpigap,pfpigap_pr,nodes)
    call mxCopyInteger1ToPtr(converged,converged_pr,1)
    call mxCopyInteger1ToPtr(reason,reason_pr,1)
    call mxCopyInteger1ToPtr(it,it_pr,1)
    call mxCopyReal8toPtr(dist_last,dist_last_pr,1)
    
    ! Deallocate arrays
    deallocate(Gg_grid)
    deallocate(Gs_grid)
    deallocate(Gmp_grid)
    deallocate(Gin_grid)
    deallocate(Gc_grid)
    deallocate(Ge_weight)  
    deallocate(Ge_nodes) 
    deallocate(Gu_weight)  
    deallocate(Gu_nodes) 
    deallocate(Gv_weight)  
    deallocate(Gv_nodes)     
    deallocate(pfc)
    deallocate(pfpigap)
    deallocate(Gg_gr)
    deallocate(Gmp_gr)
    deallocate(Gin_gr)
    deallocate(Gc_gr)
    
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
    mwSize i1,i2,i3,i4,i5,j1,j2,j3
    real(dp), allocatable, dimension(:,:,:,:) :: e_weightArr4,u_weightArr4,v_weightArr4
    real(dp), allocatable, dimension(:,:,:,:,:) :: pf_c_up,pf_pigap_up
    real(dp), allocatable, dimension(:,:,:,:,:) :: dist_c,dist_pigap
    real(dp) dist_max
    real(dp), dimension(npf,1) :: start,argzero
	real(dp), allocatable, dimension(:,:,:,:,:) :: yg,inp,izlb
	real(dp) :: locs,perbind
      
    ! Allocate memory    
    allocate(e_weightArr4(Oe_pts,Ou_pts,Ov_pts,Og_pts))
    allocate(u_weightArr4(Oe_pts,Ou_pts,Ov_pts,Os_pts))
    allocate(v_weightArr4(Oe_pts,Ou_pts,Ov_pts,Omp_pts))
	allocate(pf_c_up(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
	allocate(pf_pigap_up(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
	allocate(dist_c(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))  
	allocate(dist_pigap(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts)) 
	allocate(yg(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
	allocate(inp(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
	allocate(izlb(Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts))
	allocate(gpArr3(Oe_pts,Ou_pts,Ov_pts))
    
	! Calculate processes and weights
	do j1 = 1,Oe_pts
	  gpArr3(j1,:,:) = Ge_nodes(j1);
    end do
    
	do i1 = 1,Og_pts
	  do j1 = 1,Oe_pts
		e_weightArr4(j1,:,:,i1) = Ge_weight(i1,j1);
	  end do
    end do
	
	do i2 = 1,Os_pts
	  do j2 = 1,Ou_pts
		u_weightArr4(:,j2,:,i2) = Gu_weight(i2,j2);
	  end do 
    end do
	
	do i3 = 1,Omp_pts
      do j3 = 1,Ov_pts
        v_weightArr4(:,:,j3,i3) = Gv_weight(i3,j3);
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
        !$omp do collapse(5)
		do i5 = 1,Oc_pts
          do i4 = 1,Oin_pts
		    do i3 = 1,Omp_pts
		      do i2 = 1,Os_pts
                do i1 = 1,Og_pts  
			      ! State variables
			      start(1,1) = pfc(i1,i2,i3,i4,i5)
			      start(2,1) = pfpigap(i1,i2,i3,i4,i5)
                  state(1,1) = Gg_grid(i1)
                  state(2,1) = Gs_grid(i2)
                  state(3,1) = Gmp_grid(i3)
				  state(4,1) = Gin_grid(i4)
                  state(5,1) = Gc_grid(i5)
                  weightArr3 = e_weightArr4(:,:,:,i1)*u_weightArr4(:,:,:,i2)*v_weightArr4(:,:,:,i3)
                  
                  ! Find optimal policy functions on each node
			      argzero = start
			      call csolve(argzero)
				
			      ! Store updated policy functions            
			      pf_c_up(i1,i2,i3,i4,i5) = argzero(1,1)
			      pf_pigap_up(i1,i2,i3,i4,i5) = argzero(2,1)
                end do
              end do  
            end do
          end do
        end do
        !$omp end do
		call deallocateprivate
        !$omp end parallel

        ! Policy function distances
        dist_c = abs(pf_c_up - pfc);
        dist_pigap = abs(pf_pigap_up - pfpigap);
			
        ! Maximum distance
		dist_max = max(maxval(dist_c),maxval(dist_pigap));
			
        ! Update policy functions
        pfc = pf_c_up;
        pfpigap = pf_pigap_up;

        ! Find where ZLB binds
        !   ARC (1) and Output growth gap (3)
        yg = Gg_gr*pfc/(Pg*Gc_gr);
        !   Interest rate rule
        inp = Gin_gr**Prhoi*(Si*pfpigap**Pphipi*yg**Pphiy)**(1-Prhoi)*exp(Gmp_gr);
		!   Nodes where ZLB binds
		izlb = 0.d0
		where (inp <= 1.0_dp) izlb = 1.0_dp;
		!   Percent nodes binding
		perbind = 100*sum(izlb)/(Og_pts*Os_pts*Omp_pts*Oin_pts*Oc_pts);        
        
		! Stopping reasons
		if (dist_max > 0.5) then
			reason = 1;
		else if (all(pf_c_up < 0) .OR. any(pf_pigap_up < 0.5)) then
			reason = 2
		end if
		
        ! Check convergence criterion
        if (dist_max < Ptol) then 
            converged = 1;
		else if (reason > 0) then
			converged = 0;
		end if
        
        ! Iteration Information
        if (mod(it,25) == 1 .OR. converged >= 0) then
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
	deallocate(pf_c_up)
	deallocate(pf_pigap_up)
	deallocate(dist_c)
	deallocate(dist_pigap)
    deallocate(inp)
    deallocate(izlb)
	deallocate(gpArr3)
    
end subroutine script

! Nonlinear solver (based on Chris Sims csolve.m)
subroutine csolve(x)
    
    use structsolution
	use rnnoa_int
    use lin_sol_gen_int	            
    implicit none
    
    ! Inputs
    real(dp), dimension(npf,1), intent(inout) :: x
    
    ! Internal
    real(dp), dimension(npf) :: xtemp
    complex(dp), dimension(npf) :: fxtemp
    complex(dp), dimension(npf,1) :: fx,fmin
    real(dp), dimension(npf,1) :: fx2,xmin,randn1,dx0,dx
	real(dp), dimension(npf) :: randn0
    real(dp), parameter :: delta = 1e-6, alpha = 1e-3, crit = 1e-4
    mwSize, parameter :: itmax = 10
    real(dp), dimension(npf,npf) :: tvec,guessh,gradinv,zeros
    complex(dp), dimension(npf,npf) :: fguessh,grad
    real(dp) :: af0,af00,afmin,af
    mwSize j,itct,shrink,done,subDone
    real(dp) :: lambda,lambdamin,factor,dxSize
	mwSize :: iseed,randomize,i
    real(dp), dimension(npf,0) :: b0,x0
	real(dp) :: det0(2)
	
	iseed = 123456
	call rnset(iseed)
	zeros = 0.0d0;
	randomize = 0;
    tvec = 0.0_dp;
    do i=1,npf
        tvec(i,i) = delta;
    end do
    done = 0;
    call eqm(x,1,fx);
    af0 = sum(abs(fx));
    af00 = af0;
    itct = 0;
    do while (done == 0)
        xtemp(1) = x(1,1)
        xtemp(2) = x(2,1)
        guessh = spread(xtemp,2,npf)+tvec
		if (itct>3 .AND. af00-af0<crit*max(1.d0,af0) .AND. mod(itct,2)==1) then
			randomize = 1;
		else
			call eqm(guessh,npf,fguessh);
			fxtemp(1) = fx(1,1)
			fxtemp(2) = fx(2,1)
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
			dx0=sqrt(sum(x**2))/randn1;
	    end if
        lambda = 1;
        lambdamin = 1;
        fmin = fx;
        xmin = x;
        afmin = af0;
        dxSize = sqrt(sum(dx0**2))
        factor = 0.6_dp;
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
                factor = factor**0.6_dp;
                shrink = 1;
             end if
             if (abs(lambda*(1-factor))*dxSize > .1*delta) then
                lambda = factor*lambda;
             else if ((lambda > 0) .AND. (factor==0.6_dp)) then
                lambda = -0.3_dp;
             else
                subDone = 1;
             end if
          else if ((lambda >0) .AND. (af-af0 > (1-alpha)*lambda*af0)) then
             if (shrink == 1) then
                factor = factor**0.6_dp;
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

end subroutine csolve

subroutine eqm(x,ncol,Res)
   
	use structsolution
	use csolvepass
    implicit none
	
    ! Inputs
    mwSize, intent(in) :: ncol
    real(dp), dimension(npf,ncol), intent(in) :: x
    
    ! Output
    complex(dp), dimension(npf,ncol), intent(out) :: Res
    
    ! Internal
    real(dp) :: g,s,mp,in,c,cp,pigap,yf,yg,inp,i,lam,w
	real(dp), dimension(Oe_pts,Ou_pts,Ov_pts) :: cppArr3,pigappArr3,yfpArr3,lampArr3,sdfArr3
    real(dp) :: Ebond,Epc
    mwSize :: icol

	! State Values	
    g     = state(1,1);	  !Growth state current period
    s     = state(2,1);	  !Preference shock state current period
	mp    = state(3,1);   !Monetary policy shock current state
	in    = state(4,1);   !Notional interest rate last period
    c     = state(5,1);	  !Consumption last period
    
	do icol = 1,ncol  
		! Policy Function Guesses
		cp = x(1,icol);   	!Consumption policy current period  
        pigap = x(2,icol); 	!Inflation policy current period  
		!----------------------------------------------------------------------
		! Solve for variables
		!----------------------------------------------------------------------
		! ARC(1) and Output definition (2)
		yf = cp/(1-Pvarphip*(pigap-1)**2/2);
		! ARC (1) and Output growth gap (3)
		yg = g*cp/(Pg*c);
		! Nominal interest rate (4)
		inp = in**Prhoi*(Si*pigap**Pphipi*yg**Pphiy)**(1-Prhoi)*exp(mp);
		! Notional interest rate (5)
		i = max(1.0_dp,inp);
	!     i = abs(inp);
		! Inverse MUC (6)
		lam = cp-Ph*c/g;
		! Production function (7) and HH FOC Labor (8)
		w = Schi*yf**Peta*lam;
		!----------------------------------------------------------------------
		! Linear interpolation of the policy variables 
		!----------------------------------------------------------------------
        call allterp523( & 
			Og_pts,Os_pts,Omp_pts,Oin_pts,Oc_pts, &
			Gin_grid,Gc_grid, &
			inp,cp, &
			pfc,pfpigap, & 
			cppArr3,pigappArr3)
		!----------------------------------------------------------------------      
        ! Solve for variables inside expectations
		!----------------------------------------------------------------------  
		! ARC(1) and Output definition (6)
		yfpArr3 = cppArr3/(1-Pvarphip*(pigappArr3-1)**2/2);
		! Inverse MUC (5)
		lampArr3 = cppArr3-Ph*cp/gpArr3;
		! Stochastic discount factor
		sdfArr3 = Pbeta*lam/lampArr3;
		!----------------------------------------------------------------------
		! Numerical integration
		!---------------------------------------------------------------------- 
		Ebond = sum(weightArr3*sdfArr3/(gpArr3*(Ppi*pigappArr3)));
		Epc = sum(weightArr3*sdfArr3*(pigappArr3-1)*pigappArr3*(yfpArr3/yf));
		!----------------------------------------------------------------------
		! First-order conditions
		!----------------------------------------------------------------------
		! FOC bond
		Res(1,icol) = 1 - s*i*Ebond;
		! Firm Pricing
		Res(2,icol) = Pvarphip*(pigap-1)*pigap-(1-Pthetap)-Pthetap*w-Pvarphip*Epc;
    end do

end subroutine eqm  
 
subroutine allterp523(nx1,nx2,nx3,nx4,nx5, &
                      x4,x5, &
                      x4i,x5i, &
                      pf1,pf2, &
                      o1,o2)
    
    use structsolution, only : dp
    implicit none
	! Inputs
	mwSize, intent(in) ::  nx1,nx2,nx3,nx4,nx5
    real(dp), intent(in) :: x4i,x5i,x4(nx4),x5(nx5)
    real(dp), dimension(nx1,nx2,nx3,nx4,nx5), intent(in) :: pf1,pf2
	
	! Outputs
    real(dp), dimension(nx1,nx2,nx3), intent(out) :: o1,o2
    
	! Internal
    real(dp) :: s4, s5
    real(dp) :: x4i_min, x5i_min
    mwSize loc4, loc5
    real(dp), dimension(2) :: xi, xi_left, xi_right, w_2, w_1
    real(dp) :: w11,w12,w21,w22

    ! Grid intervals
    s4 = x4(2) - x4(1)
	s5 = x5(2) - x5(1)
    
	! Endogenous state variables   
    x4i_min = x4i - x4(1)
    loc4 = min(nx4-1,max(1,floor(x4i_min/s4) + 1));
	
	x5i_min = x5i - x5(1)
	loc5 = min(nx5-1,max(1,floor(x5i_min/s5) + 1));      

	xi = [x4i, x5i]
	xi_left = [x4(loc4), x5(loc5)]
	xi_right = [x4(loc4+1), x5(loc5+1)]

	w_2 = (xi - xi_left)/(xi_right - xi_left)
	w_1 = 1 - w_2
    
    w11 = w_1(1)*w_1(2)
    w12 = w_1(1)*w_2(2)
    w21 = w_2(1)*w_1(2)
    w22 = w_2(1)*w_2(2)   

	o1 =  	w11*pf1(:,:,:,loc4,loc5) &
		  + w12*pf1(:,:,:,loc4,loc5+1) &
		  + w21*pf1(:,:,:,loc4+1,loc5) &
		  + w22*pf1(:,:,:,loc4+1,loc5+1)
	o2 =  	w11*pf2(:,:,:,loc4,loc5) &
		  + w12*pf2(:,:,:,loc4,loc5+1) &
		  + w21*pf2(:,:,:,loc4+1,loc5) &
		  + w22*pf2(:,:,:,loc4+1,loc5+1)
	 
end subroutine allterp523
    
! Calculates algorithm and iteration duration    
subroutine itinfo(it,dist_max,perbind)
    
    use structsolution, only : dp
    implicit none    
    ! Inputs
    mwSize, intent(in) :: it
    real(dp), intent(in) :: dist_max,perbind
    
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
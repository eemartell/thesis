#include "fintrf.h"
#include "link_fnl_shared.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    ! Declarations
	implicit none

    ! mexFunction argument
    mwPointer plhs(*), prhs(*)    
    integer*4 nlhs, nrhs      
    
    ! Function declarations
    mwSize mxGetN, mxGetM, mxCreateDoubleMatrix
    mwpointer mxGetPr, mxCreateNumericArray, mxGetDimensions 
    double precision mxGetScalar  
    integer*4 mxClassIDFromClassName  
    
    ! Pointers to input/output mxArrays
    mwpointer Gk_grid_pr,Gz_grid_pr,Gbeta_grid_pr
    mwpointer Ge_weight_pr,Ge_nodes_pr
    mwpointer Gu_weight_pr,Gu_nodes_pr
    mwpointer pf_n_pr,pf_pi_pr,pf_i_pr
        
    ! Inputs
	mwSize :: nstate,npf,nshock,zlbflag
	double precision :: Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi
	double precision :: Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta
	double precision :: Schi,Sr,Sy,Ptol
    double precision, allocatable, dimension(:) :: Gk_grid,Gz_grid,Gbeta_grid
    double precision, allocatable, dimension(:) :: Ge_weight,Ge_nodes
    double precision, allocatable, dimension(:) :: Gu_weight,Gu_nodes
	double precision, allocatable, dimension(:,:,:) :: pf_n,pf_pi,pf_i
    
    ! Internal
    mwSize :: nodes
    mwSize :: Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts
    integer*4 myclassid,myclassid2,myclassid3
    
    ! Output
    mwpointer :: pf1_pr,pf2_pr,pf3_pr,converged_pr,reason_pr,it_last_pr,dist_last_pr
    double precision, allocatable, dimension(:,:,:) :: pf1,pf2,pf3
	mwSize :: converged,reason,it_last
	double precision :: dist_last

    ! Load Inputs
    ! Parameters
    nstate = mxGetScalar(prhs(1))
    npf = mxGetScalar(prhs(2))
	nshock = mxGetScalar(prhs(3))
	zlbflag = mxGetScalar(prhs(4))
    
    ! Parameters
	Palpha = mxGetScalar(prhs(5))
	Ppi = mxGetScalar(prhs(6))
	Pphipi = mxGetScalar(prhs(7))
	Pphiy = mxGetScalar(prhs(8))
	Psigma = mxGetScalar(prhs(9))
	Peta = mxGetScalar(prhs(10))
	Pvarphi = mxGetScalar(prhs(11))
	Pdelta = mxGetScalar(prhs(12))
	Pnu = mxGetScalar(prhs(13))
	Pzbar = mxGetScalar(prhs(14))
	Prhoz = mxGetScalar(prhs(15))
	Pbeta = mxGetScalar(prhs(16))
    Prhobeta = mxGetScalar(prhs(17))
	Ptheta = mxGetScalar(prhs(18))
	Schi = mxGetScalar(prhs(19))
	Sr = mxGetScalar(prhs(20))
    Sy = mxGetScalar(prhs(21))
	Ptol = mxGetScalar(prhs(22))
	
    ! Grids
    Ok_pts = mxGetN(prhs(23))
    allocate(Gk_grid(Ok_pts))
    Gk_grid_pr = mxGetPr(prhs(23))
    call mxCopyPtrToReal8(Gk_grid_pr,Gk_grid,Ok_pts)
	
    Oz_pts = mxGetN(prhs(24))
    allocate(Gz_grid(Oz_pts))
    Gz_grid_pr = mxGetPr(prhs(24))
    call mxCopyPtrToReal8(Gz_grid_pr,Gz_grid,Oz_pts)
    
    Obeta_pts = mxGetN(prhs(25))
    allocate(Gbeta_grid(Obeta_pts))
    Gbeta_grid_pr = mxGetPr(prhs(25))
    call mxCopyPtrToReal8(Gbeta_grid_pr,Gbeta_grid,Obeta_pts)         
    
    Oe_pts = mxGetM(prhs(26))
	allocate(Ge_weight(Oe_pts))
    Ge_weight_pr = mxGetPr(prhs(26))
    call mxCopyPtrToReal8(Ge_weight_pr,Ge_weight,Oe_pts)
    
    allocate(Ge_nodes(Oe_pts))
    Ge_nodes_pr = mxGetPr(prhs(27))
    call mxCopyPtrToReal8(Ge_nodes_pr,Ge_nodes,Oe_pts)	    

    Ou_pts = mxGetM(prhs(28))
	allocate(Gu_weight(Ou_pts))
    Gu_weight_pr = mxGetPr(prhs(28))
    call mxCopyPtrToReal8(Gu_weight_pr,Gu_weight,Ou_pts)
    
    allocate(Gu_nodes(Ou_pts))
    Gu_nodes_pr = mxGetPr(prhs(29))
    call mxCopyPtrToReal8(Gu_nodes_pr,Gu_nodes,Ou_pts)       
    
    ! Policy functions
    nodes = Ok_pts*Oz_pts*Obeta_pts
    allocate(pf_n(Ok_pts,Oz_pts,Obeta_pts))
    allocate(pf_pi(Ok_pts,Oz_pts,Obeta_pts))
    allocate(pf_i(Ok_pts,Oz_pts,Obeta_pts))  
    pf_n_pr = mxGetPr(prhs(30))
    pf_pi_pr = mxGetPr(prhs(31))
    pf_i_pr = mxGetPr(prhs(32))    
    call mxCopyPtrToReal8(pf_n_pr,pf_n,nodes) 
    call mxCopyPtrToReal8(pf_pi_pr,pf_pi,nodes)
    call mxCopyPtrToReal8(pf_i_pr,pf_i,nodes)            
    
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    myclassid2 = mxClassIDFromClassName('int8')
    myclassid3 = mxClassIDFromClassName('int32')    
    allocate(pf1(Ok_pts,Oz_pts,Obeta_pts))
    allocate(pf2(Ok_pts,Oz_pts,Obeta_pts))
    allocate(pf3(Ok_pts,Oz_pts,Obeta_pts))
    plhs(1) = mxCreateNumericArray(nstate,[Ok_pts,Oz_pts,Obeta_pts],myclassid,0)
    plhs(2) = mxCreateNumericArray(nstate,[Ok_pts,Oz_pts,Obeta_pts],myclassid,0)
    plhs(3) = mxCreateNumericArray(nstate,[Ok_pts,Oz_pts,Obeta_pts],myclassid,0)
	plhs(4) = mxCreateNumericArray(1,1,myclassid2,0)
	plhs(5) = mxCreateNumericArray(1,1,myclassid2,0)
	plhs(6) = mxCreateNumericArray(1,1,myclassid3,0)
	plhs(7) = mxCreateNumericArray(1,1,myclassid,0)
    pf1_pr = mxGetPr(plhs(1))
    pf2_pr = mxGetPr(plhs(2))
    pf3_pr = mxGetPr(plhs(3))
    converged_pr = mxGetPr(plhs(4))
    reason_pr = mxGetPr(plhs(5))
    it_last_pr = mxGetPr(plhs(6))
    dist_last_pr = mxGetPr(plhs(7))

    ! Run script
    call script(pf1,pf2,pf3, &
				nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts, &
				Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi, &
				Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta, &
				Schi,Sr,Sy,Ptol, &
				Gk_grid,Gz_grid,Gbeta_grid, &
                Ge_weight,Ge_nodes,Gu_weight,Gu_nodes, &
				pf_n,pf_pi,pf_i,converged,reason,it_last,dist_last)
    
    ! Load Fortran array to pointer (output to MATLAB)   
    call mxCopyReal8toPtr(pf1,pf1_pr,nodes)
    call mxCopyReal8toPtr(pf2,pf2_pr,nodes)
    call mxCopyReal8toPtr(pf3,pf3_pr,nodes)
    call mxCopyInteger1ToPtr(converged,converged_pr,1)
    call mxCopyInteger1ToPtr(reason,reason_pr,1)
    call mxCopyInteger1ToPtr(it_last,it_last_pr,1)
    call mxCopyReal8toPtr(dist_last,dist_last_pr,1)
    
    ! Deallocate arrays
    deallocate(Gk_grid)
    deallocate(Gz_grid)
    deallocate(Gbeta_grid)
    deallocate(Ge_weight)
    deallocate(Ge_nodes)
    deallocate(Gu_weight)
    deallocate(Gu_nodes)   
    deallocate(pf_n)
    deallocate(pf_pi)
    deallocate(pf_i)
    deallocate(pf1)
    deallocate(pf2)
    deallocate(pf3)
    
end subroutine mexFunction

subroutine script(pf1,pf2,pf3, &
				nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts, &
				Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi, &
				Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta, &
				Schi,Sr,Sy,Ptol, &
				Gk_grid,Gz_grid,Gbeta_grid, &
                Ge_weight,Ge_nodes,Gu_weight,Gu_nodes, &
				pf_n,pf_pi,pf_i,converged,reason,it_last,dist_last)
        
    ! Declarations
    implicit none

    ! Inputs
    mwSize :: nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts
    double precision :: Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi
	double precision :: Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta
    double precision :: Schi,Sr,Sy,Ptol
    double precision, dimension(Ok_pts) :: Gk_grid
    double precision, dimension(Oz_pts) :: Gz_grid
    double precision, dimension(Obeta_pts) :: Gbeta_grid
    double precision, dimension(Oe_pts) :: Ge_weight,Ge_nodes
    double precision, dimension(Ou_pts) :: Gu_weight,Gu_nodes
    double precision, dimension(Ok_pts,Oz_pts,Obeta_pts) :: pf_n,pf_pi,pf_i
    
    ! Internal
    mwSize i1,i2,i3,j,it
    double precision, dimension(Ok_pts,Oz_pts,Obeta_pts) :: pf_n_up,pf_pi_up,pf_i_up
    double precision, dimension(Ok_pts,Oz_pts,Obeta_pts) :: dist_n,dist_pi,dist_i
    double precision dist_max(500000)
    double precision, parameter :: tol = 1e-4, h = 1e-6
    integer, parameter :: it_max = 10
    double precision, dimension(nstate,1) :: state
    double precision, dimension(npf,1) :: start,argzero
	double precision :: diff1(50), diff2(49)
    
    ! Output
    double precision, dimension(Ok_pts,Oz_pts,Obeta_pts) :: pf1,pf2,pf3
	mwSize converged,reason,it_last
	double precision :: dist_last
	  
	it = 1;                                 ! Iteration Counter
    converged = -1;                         ! Convergence Flag
	reason = 0; 							! Stopping reason
    ! Time iteration/linear interpolation algorithm
    do while (converged == -1)
        !$omp parallel default(shared) private(start,state,argzero) 
		!$omp do collapse(3)
        do i3 = 1,Obeta_pts
		    do i2 = 1,Oz_pts
			    do i1 = 1,Ok_pts
				    ! State variables
				    start(1,1) = pf_n(i1,i2,i3)
				    start(2,1) = pf_pi(i1,i2,i3)
				    start(3,1) = pf_i(i1,i2,i3)
				    state(1,1) = Gk_grid(i1)
				    state(2,1) = Gz_grid(i2)
                    state(3,1) = Gbeta_grid(i3)
				    ! Find optimal policy functions on each node.
				    ! csolve finds the zeros of 'eqm'
				    ! Start csolve with the current policy function 
				    argzero = 0.d0
				    call csolve(start,state, & 
					    nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts, &
					    Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi, &
					    Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta, &
					    Schi,Sr,Sy, &
					    Gk_grid,Gz_grid,Gbeta_grid, &
                        Ge_weight,Ge_nodes,Gu_weight,Gu_nodes, &
					    pf_n,pf_pi,pf_i,argzero)
				 
				    ! Store updated policy functions            
				    pf_n_up(i1,i2,i3) = argzero(1,1)
				    pf_pi_up(i1,i2,i3) = argzero(2,1)
				    pf_i_up(i1,i2,i3) = argzero(3,1)
                end do
			end do
		end do
        !$omp end do
        !$omp end parallel

        ! Policy function distances    
        dist_n = abs(pf_n_up - pf_n);
        dist_pi = abs(pf_pi_up - pf_pi);
        dist_i = abs(pf_i_up - pf_i);
        
        ! Maximum distance
        dist_max(it) = max(maxval(dist_n),maxval(dist_pi),maxval(dist_i));
			
        ! Update policy functions        
        pf_n = pf_n_up;
        pf_pi = pf_pi_up;
        pf_i = pf_i_up;
                
		! Differences to check divergence
		if (it> 51) then
			diff1 = dist_max(it-50:it)-dist_max(it-51:it-1);
			diff2 = diff1(2:50)-diff1(1:49);
		end if
		
		! Stopping reasons
		if (it == 500000 &
			.OR. (it .GE. 20000 .AND. dist_max(it) > 0.001) &
			.OR. (it .GE. 100000 .AND. dist_max(it) > 0.0001)) then
			reason = 1;
		else if (all(pf_pi_up < .5) .OR. any(pf_n_up < 0)) then
			reason = 2;
		else if (it>51 .AND. all(diff1 .GE. 0) .AND. all(diff2 .GE. 0)) then
			reason = 3;
		end if
		
        ! Check convergence criterion
        if ((it > 11) .AND. all(dist_max(it-10:it) < Ptol)) then 
            converged = 1;
		else if (reason > 0) then
			converged = 0;
		end if
        
        ! Iteration Information
        if (mod(it,2) == 1 .OR. converged == 1 .OR. converged == 0) then
            call itinfo(it,dist_max(it));
        else
            it = it + 1
        end if
    end do
    
    pf1 = pf_n
    pf2 = pf_pi
    pf3 = pf_i
	it_last = it - 1;
	dist_last = dist_max(it_last);
    
end subroutine script

! Root finder based on Chris Sims csolve.m, real only, no verbose
subroutine csolve(x,state, & 
				nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts, &
				Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi, &
				Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta, &
				Schi,Sr,Sy, &
				Gk_grid,Gz_grid,Gbeta_grid, &
                Ge_weight,Ge_nodes,Gu_weight,Gu_nodes, &
				pf_n,pf_pi,pf_i,argzero)
    
    use rnnoa_int
	use rnset_int
    use cond_int
    use eye_int
    use lin_sol_gen_int
            
    implicit none
    ! Inputs
    mwSize :: nstate,npf,ncol,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts
    double precision, dimension(npf,1) :: x
    double precision, dimension(nstate,1) :: state
    double precision :: Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi
	double precision :: Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta
    double precision :: Schi,Sr,Sy
    double precision, dimension(Ok_pts) :: Gk_grid
    double precision, dimension(Oz_pts) :: Gz_grid
    double precision, dimension(Obeta_pts) :: Gbeta_grid
    double precision, dimension(Oe_pts) :: Ge_weight,Ge_nodes
    double precision, dimension(Ou_pts) :: Gu_weight,Gu_nodes
    double precision, dimension(Ok_pts,Oz_pts,Obeta_pts) :: pf_n,pf_pi,pf_i
    
    ! Internal
    double precision, dimension(npf) :: xtemp
    complex*16, dimension(npf) :: fxtemp
    complex*16, dimension(npf,1) :: fx,fmin
    double precision, dimension(npf,1) :: fx2,xmin,randn1,dx0,dx
	double precision, dimension(npf) :: randn0
    double precision, parameter :: delta = 1e-6, alpha = 1e-3, crit = 1e-4
    mwSize, parameter :: itmax = 10
    double precision, dimension(npf,npf) :: tvec,guessh,gradinv,zeros
    complex*16, dimension(npf,npf) :: fguessh,grad
    double precision :: af0,af00,afmin,af
    mwSize j,itct,shrink,done,subDone
    double precision :: lambda,lambdamin,factor,dxSize
	mwSize :: iseed,randomize
    double precision, dimension(npf,0)    :: b0,x0
	double precision :: det0(2)
    
    ! Output
    double precision, dimension(npf,1) :: argzero
	
	iseed = 123456
	call rnset(iseed)
	zeros = 0.0d0;
	randomize = 0;
    tvec = delta*eye(npf);
    done = 0;
    call eqm(x,state,1, &
			nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts, &
			Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi, &
			Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta, &
			Schi,Sr,Sy, &
			Gk_grid,Gz_grid,Gbeta_grid, &
            Ge_weight,Ge_nodes,Gu_weight,Gu_nodes, &
			pf_n,pf_pi,pf_i,fx);
    af0 = sum(abs(fx));
    af00 = af0;
    itct = 0;
    do while (done == 0)
        xtemp(1) = x(1,1)
        xtemp(2) = x(2,1)
        xtemp(3) = x(3,1)
        guessh = spread(xtemp,2,npf)+tvec
		if (itct>3 .AND. af00-af0<crit*max(1.d0,af0) .AND. mod(itct,2)==1) then
			randomize = 1;
		else
			call eqm(guessh,state,npf, &
					nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts, &
					Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi, &
					Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta, &
					Schi,Sr,Sy, &
					Gk_grid,Gz_grid,Gbeta_grid, &
                    Ge_weight,Ge_nodes,Gu_weight,Gu_nodes, &
					pf_n,pf_pi,pf_i,fguessh);
			fxtemp(1) = fx(1,1)
			fxtemp(2) = fx(2,1)
			fxtemp(3) = fx(3,1)
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
            call eqm(x+dx,state,1, &
				nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts, &
				Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi, &
				Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta, &
				Schi,Sr,Sy, &
				Gk_grid,Gz_grid,Gbeta_grid, &
                Ge_weight,Ge_nodes,Gu_weight,Gu_nodes, &
				pf_n,pf_pi,pf_i,fx);
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

subroutine eqm( x,state,ncol, &
				nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts, &
				Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi, &
				Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta, &
				Schi,Sr,Sy, &
				Gk_grid,Gz_grid,Gbeta_grid, &
                Ge_weight,Ge_nodes,Gu_weight,Gu_nodes, &
				pf_n,pf_pi,pf_i,Res)
   
    implicit none
    ! Inputs
    mwSize :: ncol,nstate,npf,nshock,zlbflag,Ok_pts,Oz_pts,Obeta_pts,Oe_pts,Ou_pts
    double precision, dimension(npf,ncol) :: x
    double precision, dimension(nstate,1) :: state
    double precision :: Palpha,Ppi,Pphipi,Pphiy,Psigma,Peta,Pvarphi
	double precision :: Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta
    double precision :: Schi,Sr,Sy
    double precision, dimension(Ok_pts) :: Gk_grid
    double precision, dimension(Oz_pts) :: Gz_grid
    double precision, dimension(Obeta_pts) :: Gbeta_grid
    double precision, dimension(Oe_pts) :: Ge_weight,Ge_nodes
    double precision, dimension(Ou_pts) :: Gu_weight,Gu_nodes
    double precision, dimension(Ok_pts,Oz_pts,Obeta_pts) :: pf_n,pf_pi,pf_i
    
    ! Output
    complex*16, dimension(npf,ncol) :: Res
    
    ! Internal
    double precision :: k,z,beta,n,pie,i
	complex*16 :: y,r,q,ytil,kac,c,w,psi,kp
    double precision, dimension(Oe_pts) :: zpVec
    double precision, dimension(Ou_pts) :: betapVec
    double precision, dimension(Oe_pts,Ou_pts) :: npMat,pipMat,ipMat
	complex*16, dimension(Oe_pts,Ou_pts) :: ypMat,qpMat,ytilpMat,kacpMat,cpMat,wpMat,psipMat,rkpMat,sdfMat
    double precision, dimension(Oe_pts,Ou_pts) :: zpMat,betapMat,eMat
    complex*16, dimension(Ou_pts) :: EfpVec,EbondVec,EconsVec
    complex*16 :: Efp,Ebond,Econs
    mwSize :: j,i1,i2
    double precision :: fastrap
	
    ! Declare local constant Pi
    integer, parameter :: DoubleReal_K = selected_real_kind (15)
    real (DoubleReal_K) :: sqrtpi
	
    ! Specify Parameters
    sqrtpi = 1.772453850905516_DoubleReal_K   

	! State Values
	k = state(1,1);     !Capital state last period  
	z = state(2,1);     !Technology state current period
    beta = state(3,1);  !Discount factor state current period
    
	do j = 1,ncol   
		! Policy Function Guesses
		n = x(1,j);     !Labor policy current period 
		pie = x(2,j);   !Inflation policy current period 
		i = x(3,j); 	!Investment cost policy current period 
		!----------------------------------------------------------------------
		! Solve for variables
		!----------------------------------------------------------------------
		! Production function
		y = z*k**Palpha*n**(1-Palpha);
		! Interest rate rule
		r = Sr*(pie/Ppi)**Pphipi*(y/Sy)**Pphiy;
		if (zlbflag == 1) then
			r = max(1.d0,abs(r));
		end if
		! Tobin's q
		q = 1+Pnu*(i/k-Pdelta);    
		! Aggregate resource constraint
		ytil = y*(1-(Pvarphi*(pie/Ppi-1)**2)/2);
		kac = Pnu*(i/k-Pdelta)**2/2;
		c = ytil - i - kac*k;    
		! FOC Labor
		w = Schi*n**Peta*c;
		! Firm FOC labor
		psi = w*n/((1-Palpha)*y);
		! Investment
		kp = i + (1-Pdelta)*k;         
		!----------------------------------------------------------------------
		! Linear interpolation of the policy variables 
		!---------------------------------------------------------------------- 
		do i2 = 1,Ou_pts
  		  do i1 = 1,Oe_pts
			zpVec(i1) = Pzbar*(z/Pzbar)**Prhoz*exp(Ge_nodes(i1))
			betapVec(i2) = Pbeta*(beta/Pbeta)**Prhobeta*exp(Gu_nodes(i2))
			npMat(i1,i2) = fastrap(Ok_pts,Oz_pts,Obeta_pts,Gk_grid,Gz_grid,Gbeta_grid,kp,zpVec(i1),betapVec(i2),pf_n)
			pipMat(i1,i2) = fastrap(Ok_pts,Oz_pts,Obeta_pts,Gk_grid,Gz_grid,Gbeta_grid,kp,zpVec(i1),betapVec(i2),pf_pi)
            ipMat(i1,i2) = fastrap(Ok_pts,Oz_pts,Obeta_pts,Gk_grid,Gz_grid,Gbeta_grid,kp,zpVec(i1),betapVec(i2),pf_i)
		  end do
		end do 
		!----------------------------------------------------------------------        
		! Solve for variables inside expectations
		!----------------------------------------------------------------------    
		zpMat = spread(zpVec,2,Ou_pts);
		betapMat = transpose(spread(betapVec,2,Oe_pts)); 		        
        ypMat = zpMat*(kp**Palpha)*(npMat**(1-Palpha));      
		qpMat = 1+Pnu*(ipMat/kp-Pdelta);  
		ytilpMat = ypMat*(1-(Pvarphi*(pipMat/Ppi-1)**2)/2);
		kacpMat = Pnu*(ipMat/kp-Pdelta)**2/2;
		cpMat = ytilpMat - ipMat - kacpMat*kp;      
		wpMat = Schi*npMat**Peta*cpMat;   
		psipMat = wpMat*npMat/((1-Palpha)*ypMat);
		rkpMat = Palpha*psipMat*ypMat/kp;
		sdfMat = betapMat*(c/cpMat);   
		!----------------------------------------------------------------------
		! Numerical integration
		!----------------------------------------------------------------------    
		eMat = spread(Ge_weight,2,Ou_pts)
        ! Apply GH across e
        EfpVec =  sum(eMat*(sdfMat*(pipMat/Ppi-1)*(ypMat/y)*(pipMat/Ppi)),1)/sqrtpi;
        EbondVec = sum(eMat*(r*sdfMat/pipMat),1)/sqrtpi;
		EconsVec = sum(eMat*(sdfMat*(rkpMat-kacpMat+(qpMat-1)*(ipMat/kp)+(1-Pdelta)*qpMat)),1)/sqrtpi;        
        ! Apply GH across u
        Efp = sum(Gu_weight*EfpVec)/sqrtpi;
		Ebond = sum(Gu_weight*EbondVec)/sqrtpi;
		Econs = sum(Gu_weight*EconsVec)/sqrtpi;		
		!----------------------------------------------------------------------
		! First-order conditions
		!----------------------------------------------------------------------
		! Firm Pricing Equation
		Res(1,j) = Pvarphi*(pie/Ppi-1)*pie/Ppi-(1-Ptheta)-Ptheta*psi-Pvarphi*Efp;
		! Bond Euler Equation
		Res(2,j) = 1-Ebond;
		! Consumption Euler Equation
		Res(3,j) = q-Econs;
	end do
	
	
end subroutine eqm
        
function fastrap(nx1,nx2,nx3,x1,x2,x3,x1i,x2i,x3i,z)

    implicit none
    mwSize nx1,nx2,nx3
    double precision :: x1i,x2i,x3i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx1,nx2,nx3) :: z
    
    double precision :: s1, s2, s3
    double precision :: x1i_min, x2i_min, x3i_min
    mwSize loc1, loc2, loc3
    double precision, dimension(3) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2, w3
    mwSize m1, m2, m3
    
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

! Calculates algorithm and iteration duration    
subroutine itinfo(it,dist_max)

    ! Inputs
    double precision :: dist_max
    
    ! Input/Output
    mwSize it
    
    ! Text Display
    integer*4, external :: mexprintf,mexEvalString
    integer*4 pm,evalout
    character*100 line
    
    ! Output iteration inforVecion
    write(line,'(A,I6,A,F17.14)') &
        'iter: ', it, &
        ' dist: ', dist_max
    pm = mexPrintf(line//achar(13)) 
    evalout = mexEvalString("pause(.0001);")

    ! Increase iteration counter
    it = it + 1;         

end subroutine itinfo
% Canonical New Keynesian Model (Rotemberg Pricing)  
%   -Habit Persistence
%   -Capital
%   -Interest Rate Smoothing
%
% Unless otherwise noted, this script and the functions it depends on are
% authored by:
%   Alexander Richter (arichter@auburn.edu)
%   Federal Reserve Bank of Dallas
%   Nathaniel Throckmorton (nathrockmorton@wm.edu)
%   College of William & Mary
%
% Copyright 2010-2016: You may freely copy and alter the following code,
% (including any scripts, functions and documentation authored by us that
% this script depends on) subject to the following conditions:
%    1) This code and any modifications to it are not sold or
%       traded in exchange for any goods or services
%    2) Credit is given to the original authors upon redistribution of
%       the original or modified code
%    3) This copyright notice is included on the original code and
%       subsequent modifications

% mex Fsolution_mod.f90 COMPFLAGS="/Qopenmp $COMPFLAGS" ...
%  -I"C:\Program Files (x86)\VNI\imsl\fnl600\intel64\include\dll"
% mex Fsolution.f90 COMPFLAGS="/Qopenmp $COMPFLAGS" ...
%     -I"C:\Program Files (x86)\VNI\imsl\fnl701\Intel64\include\dll"

clear all
% clc
tstart = tic;                           % Job timer start
%% Initialize

% Problem size
O.nstates = 6;
O.npfs = 4;
O.nshocks = 2;

% Implementation
%   F: Fortran
%   M: MATLAB
O.imp = 'M';

% Save rule {'on','off'}
saving = 'off';

% Load parameters, steady state and grids
P = parameters;
S = steadystate(P);

% Specify grid options
%   Bounds
O.cbound = [0.94,1.06]*S.c;
O.kbound = [0.85,1.25]*S.k;
O.ibound = [0.70,1.40]*S.i;
O.rnbound = [0.98,1.02]*S.r;
%   Density
% O.c_pts = 7;          % Consumption
% O.k_pts = 7;          % Capital
% O.i_pts = 7;          % Investment
% O.rn_pts = 7;         % Notional interest rate
% O.g_pts = 5;          % Growth state
% O.beta_pts = 9;       % Discount factor state
% O.mp_pts = 5;         % MP shock state
O.c_pts = 5;          % Consumption
O.k_pts = 5;          % Capital
O.i_pts = 5;          % Investment
O.rn_pts = 5;         % Notional interest rate
O.g_pts = 3;%5;          % Growth state
O.beta_pts = 5;%9;       % Discount factor state
O.mp_pts = 3;%5;         % MP shock state
O.e_pts = O.g_pts;    % Technology shock
O.u_pts = O.beta_pts; % Discount factor shock
O.v_pts = O.mp_pts;   % MP shock

% Grids
G = grids(O,P);

%% Run Policy Function Iteration Algorithm

% Obtain Guess
pf = guess(P,S,G);

% Solve the model with Fortran
if strcmp(O.imp,'F')
    disp('Solving the model with Fortran...'); pause(0.5)
    [pf.n(:),pf.q(:),pf.pi(:),pf.psi(:),R.converged,R.reason,R.it_last,R.dist_last] = ...
        Fsolution(  O.nstates,O.npfs,O.nshocks,...
                    P.tol,P.varphi,P.pi,P.zlb,P.phipi,P.phiy,...                   
                    P.theta,P.rhor,P.eta,P.h,P.alpha,P.delta,P.nu,P.g,...
                    S.r,S.chi,...
                    G.c_grid,G.k_grid,G.i_grid,G.rn_grid,...
                    G.g_grid,G.beta_grid,G.mp_grid,...
                    G.e_weight,G.e_nodes,...
                    G.u_weight,G.u_nodes,...
                    G.v_weight,G.v_nodes,...
                    pf.n,pf.q,pf.pi,pf.psi,...
                    G.c_gr,G.k_gr,G.i_gr,G.rn_gr,G.g_gr,G.mp_gr);
	toc(tstart);
% Solve the model with MATLAB
elseif strcmp(O.imp,'M')
    disp('Solving the model with MATLAB...'); pause(0.5)
    % Exogenous processes
    gpArr = G.e_nodes(:,ones(O.u_pts,1),ones(O.v_pts,1));   
    betapArr = permute(G.u_nodes(:,ones(O.e_pts,1),ones(O.v_pts,1)),[2,1,3]);
    
    % Preallocate arrays to store policy function updates
    pf_n_up = zeros(G.griddim);
    pf_q_up = zeros(G.griddim);
    pf_pi_up = zeros(G.griddim);
    pf_psi_up = zeros(G.griddim);
    it = 1;                                 % Iteration Counter
    converged = -1;                         % Convergence Flag
    reason = 0; 							% Stopping reason
    dist_max = 0;                           % Max distance vector
    while converged == -1
        istart = tic;                       % Iteration timer start
        parfor inode = 1:G.nodes
            % Find optimal policy functions on each node.
            % csolve finds the zeros of 'eqm'
            % Start csolve with the current policy function  
            start = [pf.n(inode),pf.q(inode),pf.pi(inode),pf.psi(inode)]';
            state = [G.c_gr(inode),G.k_gr(inode),G.i_gr(inode),...
                     G.rn_gr(inode),G.g_gr(inode),G.mp_gr(inode)];
            e_weightVec = G.e_weight(G.g_gr(inode) == G.g_grid,:)';
            u_weightVec = G.u_weight(G.beta_gr(inode) == G.beta_grid,:)';
            v_weightVec = G.v_weight(G.mp_gr(inode) == G.mp_grid,:)';
            e_weightArr = e_weightVec(:,ones(O.u_pts,1),ones(O.v_pts,1));   
            u_weightArr = permute(u_weightVec(:,ones(O.e_pts,1),ones(O.v_pts,1)),[2,1,3]);
            v_weightArr = permute(v_weightVec(:,ones(O.e_pts,1),ones(O.u_pts,1)),[2,3,1]);
            weightArr = e_weightArr.*u_weightArr.*v_weightArr;
            % Approximate solution
            argzero = csolve('eqm',start,[],1e-4,10,state,...
                             P,S,G,pf,...
                             gpArr,betapArr,weightArr);
            % Store updated policy functions       
            pf_n_up(inode) = argzero(1);
            pf_q_up(inode) = argzero(2);
            pf_pi_up(inode) = argzero(3);
            pf_psi_up(inode) = argzero(4);
        end

        % Policy function distances
        dist_n = abs(pf_n_up - pf.n);
        dist_q = abs(pf_q_up - pf.q);
        dist_pi = abs(pf_pi_up - pf.pi);
        dist_psi = abs(pf_psi_up - pf.psi);
        
        % Maximum distance
        dist_max = max([dist_n(:)' dist_q(:)' dist_pi(:)' dist_psi(:)']);

        % Update policy functions
        pf.n = pf_n_up;
        pf.q = pf_q_up;
        pf.pi = pf_pi_up;
        pf.psi = pf_psi_up;
        
        % Find where ZLB binds
        %   Production function
        y = (G.k_gr./G.g_gr).^P.alpha.*pf.n.^(1-P.alpha);
        % Real GDP
        rgdp = (1-(P.varphi*(pf.pi/P.pi-1).^2)/2).*y;
        rgdp_gr = G.c_gr+G.i_gr;
        %   Interest rate rule
        rn = G.rn_gr.^P.rhor.*(S.r*(pf.pi/P.pi).^P.phipi.*...
            (G.g_gr.*rgdp./(P.g*rgdp_gr)).^P.phiy).^(1-P.rhor).*exp(G.mp_gr);
        locs = find(rn <= P.zlb);
        %   Percent nodes binding
        perbind = 100*numel(locs)/G.nodes;

        % Stopping reasons
        if dist_max > 0.1
            reason = 1;
        elseif (all(pf_n_up(:) < 0) || any(pf_pi_up(:) < 0.5))
            reason = 2;
        end

        % Check convergence criterion
        if dist_max < P.tol
            converged = 1;
        elseif reason > 0
            converged = 0;
        end

        % Iteration Information
        if mod(it,1) == 0 || converged == 1
            it = itinfo(istart,tstart,1,it,dist_max);
            if perbind>0
                disp(['nodes binding: ' num2str(perbind) '%']);
            end
        else
            it = it+1;
        end
    end
else
    load('solution.mat')
end

%% Save results
if strcmp(saving,'on')
    save('solutions/solution_test.mat','pf','O','P','S','G','R');
end
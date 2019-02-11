% Time Iteration (Linear Interpolation):
% Canonical New Keynesian Model (Rotemberg Pricing) with Capital  
%   -Interest rate rule responds to inflation and output
%   -Balanced budget FP (no debt)
%   -Imposes the zero lower bound on the interest rate
%   -Stochastic technology process
%
% Unless otherwise noted, this script and the functions it depends on are
% authored by:
%   Alexander Richter (arichter@auburn.edu)
%   Auburn University - Department of Economics
%   Nathaniel Throckmorton (nathrock@indiana.edu)
%   Indiana University - Department of Economics
%
% Copyright 2010-2014: You may freely copy and alter the following code,
% (including any scripts, functions and documentation authored by us that
% this script depends on) under the following conditions:
%    1) This code and any modifications to it are not sold or
%       traded in exchange for any goods or services
%    2) Credit is given to the original authors upon redistribution of
%       the original or modified code
%    3) This copyright notice is included on the original code and
%       subsequent modifications
% mex Fscript.f90 COMPFLAGS="/Qopenmp $COMPFLAGS" ...
%  -I"C:\Program Files (x86)\VNI\imsl\fnl600\intel64\include\dll"

clear
% clc
tstart = tic;                           % Job timer start
%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------
% Problem size
O.nstates = 3;
O.npfs = 3;
O.nshocks = 2;

% ZLB flag
%  0: do not impose ZLB constraint
%  1: impose ZLB constraint
O.zlbflag = 1;

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
O.kbound = [0.95 1.05]*S.k;

pr = .000005;
sigz = sqrt(P.sige^2/(1-P.rhoz^2));
left = norminv(pr,P.zbar,sigz);
right = norminv(1-pr,P.zbar,sigz);
O.zbound = [left right];

sigbeta = sqrt(P.sigu^2/(1-P.rhobeta^2));
left = norminv(pr,P.beta,sigbeta);
right = norminv(1-pr,P.beta,sigbeta);
O.betabound = [left right];

pts = 101;
O.k_pts = pts;
O.z_pts = pts;
O.beta_pts = pts;
O.e_pts = 21;
O.u_pts = 21;

% Grids
G = grids(O,P);

% Obtain Guess
pf = guess(P,S,G);
%% ------------------------------------------------------------------------
% Solve Full Model
%--------------------------------------------------------------------------
% Apply Time Iteration Algorithm
if strcmp(O.imp,'F')
    [pf.n,pf.pi,pf.i,converged,reason,it_last,dist_last] = ...
    Fscript(O.nstates,O.npfs,O.nshocks,O.zlbflag,...
            P.alpha,P.pi,P.phipi,P.phiy,P.sigma,P.eta,P.varphi,...
            P.delta,P.nu,P.zbar,P.rhoz,P.beta,P.rhobeta,P.theta,...
            S.chi,S.r,S.y,P.tol,...
            G.k_grid,G.z_grid,G.beta_grid,...
            G.e_weight,G.e_nodes,G.u_weight,G.u_nodes, ...
            pf.n,pf.pi,pf.i);
	toc(tstart);
elseif strcmp(O.imp,'M')
    % Turn on matlabpool 
    %if matlabpool('size') == 0
    %    matlabpool
    %end
    
    % Initialize policy function updates
    pf_n_up = zeros(G.griddim);
    pf_pi_up = pf_n_up;
    pf_i_up = pf_n_up;

    tstart = tic;                           % Timer start
    it = 1;                                 % Iteration Counter
    converged = 0;                          % Convergence Flag
    dist_max = 0;                           % Max distance vector
    reason = 0;                             % non-convergence reason
    while converged == 0
        istart = tic;                       % Iteration timer start    
        for i = 1:G.nodes
            state = [G.k_gr(i) G.z_gr(i) G.beta_gr(i)];
            start = [pf.n(i) pf.pi(i) pf.i(i)]';
            argzero = csolve('eqm',start,[],1e-4,10,state, ....
                        O.zlbflag,...
                        P.alpha,P.pi,P.phipi,P.phiy,P.eta,P.varphi,...  
                        P.delta,P.nu,P.zbar,P.rhoz,P.beta,P.rhobeta,P.theta,...
                        S.chi,S.r,S.y,...       
                        G.e_weight,G.e_nodes,G.u_weight,G.u_nodes,...
                        G.k_grid,G.z_grid,G.beta_grid,...
                        pf.n,pf.pi,pf.i);        
            % Store updated rules     
            pf_n_up(i) = argzero(1);
            pf_pi_up(i) = argzero(2);
            pf_i_up(i) = argzero(3);       
        end

        % Rule distances between iterations
        dist_n = abs(pf_n_up - pf.n);
        dist_pi = abs(pf_pi_up - pf.pi);
        dist_i = abs(pf_i_up - pf.i);

        % Maximum distance across all nodes
        dist_max(it) = max([dist_n(:)' dist_pi(:)' dist_i(:)']);

        % Replace current rule with update
        pf.n = pf_n_up;
        pf.pi = pf_pi_up;
        pf.i = pf_i_up;

        % Production function
        y = G.z_gr.*G.k_gr.^P.alpha.*pf.n.^(1-P.alpha);    
        % Interest rate rule
        r = max(1,S.r*(pf.pi/P.pi).^P.phipi.*(y./S.y).^P.phiy);
        %   Percent nodes binding
        perbind = 100*sum(r(:)==1)/G.nodes;

        % Differences to check divergence
        if it > 51
            diff1 = dist_max(it-50:it)-dist_max(it-51:it-1);
            diff2 = diff1(2:50)-diff1(1:49);
        end

        % Stopping reasons
        if it == 500000 ...
            || (it >= 20000 && dist_max(it) > 0.001) ...
            || (it >= 100000 && dist_max(it) > 0.0001)
            reason = 1;
        elseif all(pf.pi(:) < .5) || any(pf.n(:) < 0)
            reason = 2;
        elseif perbind >= 50
            reason = 3;
        elseif (it>51 && all(diff1(:) >= 0) && all(diff2(:) >= 0))
            reason = 4;
        end

        % Check convergence criterion
        if it > 11 && all(dist_max(end-10:end) < P.tol)
            converged = 1;
        end

        % Iteration Information and save policy functions to disk
        if mod(it,10) == 1 || converged == 1
            it = itinfo(istart,tstart,1,it,dist_max(it));
        else
            it = it+1;
        end
    end
end

% Find where ZLB binds
% Production function
y = G.z_gr.*G.k_gr.^P.alpha.*pf.n.^(1-P.alpha); 
%   Interest rate rule
r = max(1,S.r*(pf.pi/P.pi).^P.phipi.*(y./S.y).^P.phiy);
locs = find(r == 1);
%   Percent nodes binding
perbind = 100*numel(locs)/G.nodes;
disp(['nodes binding: ' num2str(perbind) '%']);

% Save results
if strcmp(saving,'on')
    % Save baseline 
    fname = 'Rules\pf_baseline';
    save([fname '.mat'],'pf','O','P','S','G');

    % Create new results file and save baselinename
    baselinename = fname;
    save('Rules\results.mat','baselinename')
end
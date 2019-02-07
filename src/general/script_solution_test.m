% Canonical New Keynesian Model (Rotemberg Pricing)  
%
% Unless otherwise noted, this script and the functions it depends on are
% authored by:
%   Alexander W. Richter (Alex.Richter@dal.frb.org)
%   Federal Reserve Bank of Dallas
%   Nathaniel A. Throckmorton (nathrockmorton@wm.edu)
%   College of William and Mary
%
% Copyright 2010-2016, Creative Commons, BY-NC-SA 4.0
%    Attribution-NonCommercial-ShareAlike 4.0 International
%    https://creativecommons.org/licenses/by-nc-sa/4.0/

% mex Fsolution_test.f90 COMPFLAGS="/Qopenmp $COMPFLAGS" ...
%  -I"C:\Program Files (x86)\VNI\imsl\fnl701\Intel64\include\dll"
% mex Fallterp523_R.f90

clear all
% clc
tstart = tic;                           % Job timer start
%% Initialize
% Implementation
%   F: Fortran
%   M: MATLAB
imp = 'M';

% Save rule {'on','off'}
saving = 'off';

% Load options, parameters, and steady state
load('options.mat');

%% Run Policy Function Iteration Algorithm

% Obtain Guess
pf = guess(P,S,G);

% Solve the model with Fortran
if strcmp(imp,'F')
    disp('Solving the model with Fortran...'); pause(0.5)
    [pf.c(:),pf.pigap(:),R.converged,R.reason,R.it_last,R.dist_last] = Fsolution_test( ...
        O.nstates,O.npfs,O.nshocks,...
        P.tol,P.beta,P.thetap,P.eta,P.g,P.pi,P.s,...
        P.varphip,P.h,P.rhos,P.rhoi,P.phipi,P.phiy,...
        S.chi,S.i,...
        G.g_grid,G.s_grid,G.mp_grid,G.in_grid,G.c_grid,...
        G.e_weight,G.e_nodes,...
        G.u_weight,G.u_nodes,...
        G.v_weight,G.v_nodes,...
        pf.c,pf.pigap,...
        G.g_gr,G.mp_gr,G.in_gr,G.c_gr);
	R.solution_time = toc(tstart);
% Solve the model with MATLAB
elseif strcmp(imp,'M')
    disp('Solving the model with MATLAB...'); pause(0.5)
    % Exogenous processes   
    gpArr3 = repmat(G.e_nodes,[1,O.epss_pts,O.epsmp_pts]); 
    
    % Preallocate arrays to store policy function updates
    pf_c_up = zeros(G.griddim);
    pf_pigap_up = zeros(G.griddim);
    it = 1;                                 % Iteration Counter
    converged = -1;                         % Convergence Flag
    reason = 0; 							% Stopping reason
    dist_max = 0;                           % Max distance vector
    while converged == -1
        istart = tic;                       % Iteration timer start
        for inode = 1:G.nodes
            % Find optimal policy functions on each node  
            start = [pf.c(inode),pf.pigap(inode)]';
            state = [G.g_gr(inode),G.s_gr(inode),G.mp_gr(inode),G.in_gr(inode),G.c_gr(inode)]; 
            e_weightVec = G.e_weight(G.g_gr(inode) == G.g_grid,:)';
            u_weightVec = G.u_weight(G.s_gr(inode) == G.s_grid,:)';
            v_weightVec = G.v_weight(G.mp_gr(inode) == G.mp_grid,:)';
            
            e_weightArr3 = e_weightVec(:,ones(O.epss_pts,1),ones(O.epsmp_pts,1));
            u_weightArr3 = permute(u_weightVec(:,ones(O.epsg_pts,1),ones(O.epsmp_pts,1)),[2,1,3]);
            v_weightArr3 = permute(v_weightVec(:,ones(O.epsg_pts,1),ones(O.epss_pts,1)),[2,3,1]);
            weightArr3 = e_weightArr3.*u_weightArr3.*v_weightArr3;
            % Approximate solution
            argzero = csolve('eqm',start,[],1e-4,10,state,...
                             O,P,S,G,pf,gpArr3,weightArr3);
            % Store updated policy functions       
            pf_c_up(inode) = argzero(1);
            pf_pigap_up(inode) = argzero(2);
        end

        % Policy function distances
        dist_c = abs(pf_c_up - pf.c);
        dist_pigap = abs(pf_pigap_up - pf.pigap);
        
        % Maximum distance
        dist_max = max([dist_c(:)',dist_pigap(:)']);

        % Update policy functions
        pf.c = pf_c_up;
        pf.pigap = pf_pigap_up;
        
        % Find where ZLB binds
        %   ARC (1) and Output growth gap (3)
        yg = G.g_gr.*pf.c./(P.g*G.c_gr);
        %   Interest rate rule
        inp = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi.*yg.^P.phiy).^(1-P.rhoi).*exp(G.mp_gr);
        %   Percent nodes binding
        locs = find(inp <= 1);
        perbind = 100*numel(locs)/G.nodes;

        % Stopping reasons
        if dist_max > 0.5
            reason = 1;
        elseif (all(pf_c_up(:) < 0) || any(pf_pigap_up(:) < 0.5))
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
end

%% Save results
if strcmp(saving,'on')
    if P.zlbflag    
        save('solution_test.mat','pf','O','P','S','G');
    else
        save('solution_test_nozlb.mat','pf','O','P','S','G');
    end
end
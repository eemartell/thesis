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

clear all
% clc
tstart = tic;                           % Job timer start
%% Initialize

% Save rule {'on','off'}
saving = 'on';

% Load options, parameters, and steady state
load('solutions/options.mat');

% Iteration
%   ti: time iteration
%   fp: fixed point
O.it = 'fp';

%% Run Policy Function Iteration Algorithm

% Obtain Guess
pf = guess(P,S,G);

disp('Solving the model with MATLAB...'); pause(0.5)
% Exogenous processes   
gpArr3 = repmat(G.e_nodes,[1,O.u_pts,O.v_pts]); 
spArr3 = permute(repmat(G.u_nodes,[1,O.e_pts,O.v_pts]),[2,1,3]); 

% Preallocate arrays to store policy function updates
pf_c_up = zeros(G.griddim);
pf_pigap_up = zeros(G.griddim);
it = 1;                                 % Iteration Counter
converged = -1;                         % Convergence Flag
reason = 0; 							% Stopping reason
dist_max = 0;                           % Max distance vector
while converged == -1
istart = tic;                       % Iteration timer start
%        parfor inode = 1:G.nodes
for inode = 1:G.nodes
    % Find optimal policy functions on each node  
    start = [pf.c(inode),pf.pigap(inode)]';
    state = [G.g_gr(inode),G.s_gr(inode),G.mp_gr(inode),G.in_gr(inode)]; %%%no in
    e_weightVec = G.e_weight(G.g_gr(inode) == G.g_grid,:)';
    u_weightVec = G.u_weight(G.s_gr(inode) == G.s_grid,:)';
    v_weightVec = G.v_weight(G.mp_gr(inode) == G.mp_grid,:)';
    
    e_weightArr3 = e_weightVec(:,ones(O.u_pts,1),ones(O.v_pts,1));
    u_weightArr3 = permute(u_weightVec(:,ones(O.e_pts,1),ones(O.v_pts,1)),[2,1,3]);
    v_weightArr3 = permute(v_weightVec(:,ones(O.e_pts,1),ones(O.u_pts,1)),[2,3,1]);
    weightArr3 = e_weightArr3.*u_weightArr3.*v_weightArr3;
    % Approximate solution
    if strcmp(O.it,'ti')
        argzero = csolve('eqm',start,[],1e-4,10,state,...
                      O,P,S,G,pf,gpArr3,spArr3,weightArr3);
    elseif strcmp(O.it,'fp')
        argzero = eqm_fp(start,state,...
                         O,P,S,G,pf,gpArr3,spArr3,weightArr3); %%%
    end
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
inp = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr); %%%no in_gr, make like eq in new Model.pdf
locs = find(inp <= 1);
%   Percent nodes binding
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
%% Save results
if strcmp(saving,'on')
    fname = ['solution' O.it num2str(P.zlbflag)];
    save(['solutions/' fname],'pf','O','P','S','G','V','perbind');    
end

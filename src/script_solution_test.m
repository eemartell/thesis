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

% Solution algorithm
%   ART:  pf.hh = pf.hh;    pf.firm = pf.firm
%   Gust: pf.hh = pf.Vlam; pf.firm = pf.Vpi 
O.alg = 'Gust';

%% Run Policy Function Iteration Algorithm

% Obtain Guess
%%%if Gust, add pf.hh_zlb and pf.firm_zlb
pf = guess(P,S,G,O);
if strcmp(O.alg,'Gust')
    pf.hh_zlb = pf.hh;
    pf.firm_zlb = pf.firm;
end

disp('Solving the model with MATLAB...'); pause(0.5)
% Exogenous processes   
gpArr3 = repmat(G.e_nodes,[1,O.u_pts,O.v_pts]); 
spArr3 = permute(repmat(G.u_nodes,[1,O.e_pts,O.v_pts]),[2,1,3]); 
mpArr3 = permute(repmat(G.v_nodes,[1,O.e_pts,O.u_pts]),[3,1,2]); 

% Preallocate arrays to store policy function updates
%%% preallocate all 4 if Gust
pf_hh_up = zeros(G.griddim);
pf_firm_up = zeros(G.griddim);
if strcmp(O.alg,'Gust')
    pf_hh_zlb_up = zeros(G.griddim);
    pf_firm_zlb_up = zeros(G.griddim);
end
it = 1;                                 % Iteration Counter
converged = -1;                         % Convergence Flag
reason = 0; 							% Stopping reason
dist_max = 0;                           % Max distance vector
while converged == -1
istart = tic;                       % Iteration timer start
%        parfor inode = 1:G.nodes
for inode = 1:G.nodes
    % Find optimal policy functions on each node 
    if strcmp(O.alg,'ART')
        start = [pf.hh(inode),pf.firm(inode)]'; %%%unpack all 4 if Gust
    elseif strcmp(O.alg,'Gust')
        start = [pf.hh(inode),pf.hh_zlb(inode),pf.firm(inode),pf.firm_zlb(inode)]';
    end
    state = [G.g_gr(inode),G.s_gr(inode),G.mp_gr(inode),G.in_gr(inode)]; 
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
        if strcmp(O.alg,'ART')
            argzero = eqm_fp(start,state,...
                            O,P,S,G,pf,gpArr3,spArr3,weightArr3);
        elseif strcmp(O.alg,'Gust')
            argzero = eqm_fp_gustetal(start,state,...
                            O,P,S,G,pf,gpArr3,mpArr3,weightArr3);
        end
    end
    %%%store all 4 if Gust
    % Store updated policy functions
    if strcmp(O.alg,'ART')
        pf_hh_up(inode) = argzero(1);
        pf_firm_up(inode) = argzero(2);
    elseif strcmp(O.alg,'Gust')
        pf_hh_up(inode) = argzero(1);
        pf_hh_zlb_up(inode) = argzero(2);
        pf_firm_up(inode) = argzero(3);
        pf_firm_zlb_up(inode) = argzero(4);
    end
end

%%%get distances of all 4 if Gust
% Policy function distances
dist_hh = abs(pf_hh_up - pf.hh);
dist_firm = abs(pf_firm_up - pf.firm);
if strcmp(O.alg,'Gust')
    dist_hh_zlb = abs(pf_hh_zlb_up - pf.hh_zlb);
    dist_firm_zlb = abs(pf_firm_zlb_up - pf.firm_zlb);
end

%%%get max dist of all 4 if Gust
% Maximum distance
if strcmp(O.alg,'ART')
    dist_max = max([dist_hh(:)',dist_firm(:)']);
elseif strcmp(O.alg,'Gust')
    dist_max = max([dist_hh(:)',dist_hh_zlb(:)',dist_firm(:)',dist_firm_zlb(:)']);
end

%%%update all 4 if Gust
% Update policy functions
pf.hh = pf_hh_up;
pf.firm = pf_firm_up;
if strcmp(O.alg,'Gust')
    pf.hh_zlb = pf_hh_zlb_up;
    pf.firm_zlb = pf_firm_zlb_up;
end

% Find where ZLB binds
if strcmp(O.alg,'ART')
    inp = G.in_gr.^P.rhoi.*(S.i*pf_firm_up.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);
elseif strcmp(O.alg,'Gust')
    pigap_up = (1+sqrt((P.varphi + 4*pf_firm_up)/P.varphi))/2;
    inp = G.in_gr.^P.rhoi.*(S.i*pigap_up.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);
end
locs = find(inp <= 1);
%   Percent nodes binding
perbind = 100*numel(locs)/G.nodes;

% Stopping reasons
if dist_max > 0.75
    reason = 1;
end
if strcmp(O.alg,'ART')
    if (all(pf_hh_up(:) < 0) || any(pf_firm_up(:) < 0.5))
        reason = 2;
    end
elseif strcmp(O.alg,'Gust')
    c_up = 1/pf_firm_up;
    c_zlb_up = 1/pf_firm_zlb_up;
    pigap_zlb_up = (1+sqrt((P.varphi + 4*pf_firm_zlb_up)/P.varphi))/2;    
    if (all(c_up(:) < 0) || all(c_zlb_up(:) < 0) || any(pigap_up(:) < 0.5)...
        || any(pigap_zlb_up(:) < 0.5))
        reason = 2;
    end
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
    fname = ['solution' O.it O.alg];
    save(['solutions/' fname],'pf','O','P','S','G','V','perbind');    
end

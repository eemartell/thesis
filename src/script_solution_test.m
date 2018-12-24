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
pf = guess(P,S,G,O);

% Obtain guess for ZLB nodes binding
% Find where ZLB binds
if strcmp(O.alg,'ART')
    inp = G.in_gr.^P.rhoi.*(S.i*pf.firm.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);
elseif strcmp(O.alg,'Gust')
    pigap_up = (1+sqrt((P.varphi + 4*pf.firm)/P.varphi))/2;
    inp = G.in_gr.^P.rhoi.*(S.i*pigap_up.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);
end
statezlbinfo = (inp <= 1);

disp('Solving the model with MATLAB...'); pause(0.5)
% Exogenous processes   
gpArr3 = repmat(G.e_nodes,[1,O.u_pts,O.v_pts]); 
spArr3 = permute(repmat(G.u_nodes,[1,O.e_pts,O.v_pts]),[2,1,3]); 

% Preallocate arrays to store policy function updates
pf_hh_up = zeros(G.griddim);
pf_firm_up = zeros(G.griddim);
it = 1;                                 % Iteration Counter
converged = -1;                         % Convergence Flag
reason = 0; 							% Stopping reason
dist_max = 0;                           % Max distance vector
while converged == -1
istart = tic;                       % Iteration timer start
%        parfor inode = 1:G.nodes
for inode = 1:G.nodes
    % Find optimal policy functions on each node  
    start = [pf.hh(inode),pf.firm(inode)]';
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
                            O,P,S,G,pf,gpArr3,spArr3,weightArr3);
            if statezlbinfo(inode) == 1 
                % does something diff need to happen in eqm for zlb case?
                argzerozlb = eqm_fp_gustetal(start,state,...
                                O,P,S,G,pf,gpArr3,spArr3,weightArr3);
            else
                argzerozlb = argzero;
            end
        end
    end
    % Store updated policy functions       
    pf_hh_up(inode) = argzero(1);
    pf_firm_up(inode) = argzero(2);
end

% Policy function distances
dist_hh = abs(pf_hh_up - pf.hh);
dist_firm = abs(pf_firm_up - pf.firm);

% Maximum distance
dist_max = max([dist_hh(:)',dist_firm(:)']);

% Update policy functions
pf.hh = pf_hh_up;
pf.firm = pf_firm_up;

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
if dist_max > 1%0.5
    reason = 1;
end
if strcmp(O.alg,'ART')
    if (all(pf_hh_up(:) < 0) || any(pf_firm_up(:) < 0.5))
        reason = 2;
    end
elseif strcmp(O.alg,'Gust')
    c_up = 1/pf_firm_up;
    if (all(c_up(:) < 0) || any(pigap_up(:) < 0.5))
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

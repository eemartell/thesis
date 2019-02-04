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
tic
time = zeros(300,1);
%% Initialize

% Save rule {'on','off'}
saving = 'on';

% Solution algorithm
O.alg = 'ART';

% Load options, parameters, and steady state
if strcmp(O.alg,'ART')
    load('solutions/optionsART.mat');
elseif strcmp(O.alg,'Gust')
    load('solutions/optionsGust.mat');
end

% Iteration
%   ti: time iteration
%   fp: fixed point
O.it = 'fp';
%% Run Policy Function Iteration Algorithm

% Obtain Guess
pf = guess(P,S,G,O);
if strcmp(O.alg,'Gust')
    pf.c_zlb = pf.c;
end

disp('Solving the model with MATLAB...'); pause(0.5)
% Exogenous processes   
gpArr3 = repmat(G.e_nodes,[1,O.u_pts,O.v_pts]); 
apArr3 = permute(repmat(G.u_nodes,[1,O.e_pts,O.v_pts]),[2,1,3]); 
mpArr3 = permute(repmat(G.v_nodes,[1,O.e_pts,O.u_pts]),[2,3,1]); 

% Preallocate arrays to store policy function updates
pf_c_up = zeros(G.griddim);
pf_pigap_up = zeros(G.griddim);
if strcmp(O.alg,'Gust')
    pf_c_zlb_up = zeros(G.griddim);
end
it = 1;                                 % Iteration Counter
converged = -1;                         % Convergence Flag
reason = 0; 							% Stopping reason
dist_max = 0;                           % Max distance vector
toc %comment out if you don't want to see time before iterations
tic
while converged == -1
istart = tic;                       % Iteration timer start
%        parfor inode = 1:G.nodes
for inode = 1:G.nodes
    % Find optimal policy functions on each node 
    if strcmp(O.alg,'ART')
        start = [pf.c(inode),pf.pigap(inode)]';
    elseif strcmp(O.alg,'Gust')
        start = [pf.c(inode),pf.c_zlb(inode),pf.pigap(inode)]';
    end
    state = [G.g_gr(inode),G.a_gr(inode),G.mp_gr(inode),G.in_gr(inode)]; 
    e_weightVec = G.e_weight(G.g_gr(inode) == G.g_grid,:)';
    u_weightVec = G.u_weight(G.a_gr(inode) == G.a_grid,:)';
    v_weightVec = G.v_weight(G.mp_gr(inode) == G.mp_grid,:)';
    
    e_weightArr3 = e_weightVec(:,ones(O.u_pts,1),ones(O.v_pts,1));
    u_weightArr3 = permute(u_weightVec(:,ones(O.e_pts,1),ones(O.v_pts,1)),[2,1,3]);
    v_weightArr3 = permute(v_weightVec(:,ones(O.e_pts,1),ones(O.u_pts,1)),[2,3,1]);
    weightArr3 = e_weightArr3.*u_weightArr3.*v_weightArr3;
    % Approximate solution
    if strcmp(O.it,'ti')
        if strcmp(O.alg,'ART')
            argzero = csolve('eqm',start,[],1e-4,10,state,...
                      O,P,S,G,pf,gpArr3,apArr3,weightArr3);
        elseif strcmp(O.alg,'Gust')
            argzero = csolve('eqm_gustetal',start,[],1e-4,10,state,...
                      O,P,S,G,pf,gpArr3,apArr3,weightArr3);
        end

    elseif strcmp(O.it,'fp')
        if strcmp(O.alg,'ART')
            argzero = eqm_fp(start,state,...
                            O,P,S,G,pf,gpArr3,apArr3,weightArr3);
        elseif strcmp(O.alg,'Gust')
            argzero = eqm_fp_gustetal(start,state,...
                            O,P,S,G,pf,gpArr3,mpArr3,weightArr3);
        end
    end
    % Store updated policy functions
    if strcmp(O.alg,'ART')
        pf_c_up(inode) = argzero(1);
        pf_pigap_up(inode) = argzero(2);
    elseif strcmp(O.alg,'Gust')
        pf_c_up(inode) = argzero(1);
        pf_c_zlb_up(inode) = argzero(2);
        pf_pigap_up(inode) = argzero(3);
    end
end

% Policy function distances
dist_c = abs(pf_c_up - pf.c);
dist_pigap = abs(pf_pigap_up - pf.pigap);
if strcmp(O.alg,'Gust')
    dist_c_zlb = abs(pf_c_zlb_up - pf.c_zlb);
end

% Maximum distance
if strcmp(O.alg,'ART')
    dist_max = max([dist_c(:)',dist_pigap(:)']);
elseif strcmp(O.alg,'Gust')
    dist_max = max([dist_c(:)',dist_c_zlb(:)',dist_pigap(:)']);
end

% Update policy functions
pf.c = pf_c_up;
pf.pigap = pf_pigap_up;
if strcmp(O.alg,'Gust')
    pf.c_zlb = pf_c_zlb_up;
end

% Find where ZLB binds
inp = G.in_gr.^P.rhoi.*(S.i*pf_pigap_up.^P.phipi).^(1-P.rhoi).*exp(G.mp_gr);
locs = find(inp <= 1);
%   Percent nodes binding
perbind = 100*numel(locs)/G.nodes;

% Stopping reasons
if dist_max > 0.75
    reason = 1;
end
if strcmp(O.alg,'ART')
    if (all(pf_c_up(:) < 0) || any(pf_pigap_up(:) < 0.5))
        reason = 2;
    end
elseif strcmp(O.alg,'Gust')
    if (all(pf_c_up(:) < 0) || all(pf_c_zlb_up(:) < 0) || any(pf_pigap_up(:) < 0.5))
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
time(it-1) = toc;
toc %comment out if you don't want to see times on each iteration
tic
end
time = time(time>0);
%% Save results
if strcmp(saving,'on')
    fname = ['solution' O.it num2str(O.a_pts) O.alg];
    save(['solutions/' fname],'pf','O','P','S','G','V','perbind');    
end
disp(fname)

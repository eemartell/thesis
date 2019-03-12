% New Keynesian Model
%    ZLB Constraint
%    Habit persistence
%    Interest rate smoothing
%    Capital
%    Investment adjustment costs
%  	 Variable Capital Utilization
%    Wage adjustment costs
%
% Unless otherwise noted, this script and the functions it depends on are
% authored by:
%   Alexander Richter (Alex.Richter@dal.frb.org)
%   Federal Reserve Bank of Dallas
%   Nathaniel Throckmorton (nathrockmorton@wm.edu)
%   William & Mary
%
% Copyright 2010-2018: You may freely copy and alter the following code,
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

% Save rule {'on','off'}
saving = 'off';

% Load parameters, steady state and grids
load('options.mat')

% Iteration
%   ti: time iteration
%   fp: fixed point
O.it = 'fp';

% Solution algorithm
O.alg = 'ART';

%% Run Policy Function Iteration Algorithm

% Obtain Guess
pf = guess(V,P,S,G);
if strcmp(O.alg,'Gust')
    pf.n_zlb = pf.n;
end

% Exogenous processes
gpArr3 = G.epsg_nodes(:,ones(O.epss_pts,1),ones(O.epsmp_pts,1));   
spArr3 = permute(repmat(G.epss_nodes,[1,O.epsg_pts,O.epsmp_pts]),[2,1,3]); 
mpArr3 = permute(repmat(G.epsmp_nodes,[1,O.epsg_pts,O.epss_pts]),[2,3,1]); 

% Preallocate arrays to store policy function updates
pf_pigap_up = zeros(G.griddim);
pf_n_up = zeros(G.griddim);
pf_q_up = zeros(G.griddim);
pf_mc_up = zeros(G.griddim);
if strcmp(O.alg,'Gust')
    pf_n_zlb_up = zeros(G.griddim);
end

it = 1;                                 % Iteration Counter
converged = -1;                         % Convergence Flag
reason = 0; 							% Stopping reason
dist_max = 0;                           % Max distance vector
while converged == -1
    istart = tic;                       % Iteration timer start
    parfor inode = 1:G.nodes
        % Find optimal policy functions on each node
        if strcmp(O.alg,'ART')
            start = [pf.pigap(inode),pf.n(inode),pf.q(inode),pf.mc(inode)]';
        elseif strcmp(O.alg,'Gust')
            start = [pf.pigap(inode),pf.n(inode),pf.n_zlb(inode),pf.q(inode),pf.mc(inode)]';
        end
        state = [G.g_gr(inode),G.s_gr(inode),G.mp_gr(inode),...
            G.in_gr(inode),G.c_gr(inode),G.k_gr(inode),G.x_gr(inode)];
        epsg_weightVec = G.epsg_weight(G.g_gr(inode) == G.g_grid,:)';
        epss_weightVec = G.epss_weight(G.s_gr(inode) == G.s_grid,:)';
        epsmp_weightVec = G.epsmp_weight(G.mp_gr(inode) == G.mp_grid,:)';
        epsg_weightArr = epsg_weightVec(:,ones(O.epss_pts,1),ones(O.epsmp_pts,1));   
        epss_weightArr = permute(epss_weightVec(:,ones(O.epsg_pts,1),ones(O.epsmp_pts,1)),[2,1,3]);
        epsmp_weightArr = permute(epsmp_weightVec(:,ones(O.epsg_pts,1),ones(O.epss_pts,1)),[2,3,1]);
        weightArr3 = epsg_weightArr.*epss_weightArr.*epsmp_weightArr;
        % Approximate solution
        if strcmp(O.it,'ti')
            if strcmp(O.alg,'ART')
                argzero = csolve('eqm',start,[],1e-4,10,state,...
                        O,P,S,G,pf,gpArr3,weightArr3);
            elseif strcmp(O.alg,'Gust')
                argzero = csolve('eqm_gustetal',start,[],1e-4,10,state,...
                      O,P,S,G,pf,gpArr3,mpArr3,weightArr3);
            end
        elseif strcmp(O.it,'fp')
            if strcmp(O.alg,'ART')
                argzero = eqm_fp(start,state,...
                            O,P,S,G,pf,gpArr3,weightArr3);
            elseif strcmp(O.alg,'Gust')
                argzero = eqm_fp_gustetal(start,state,...
                            O,P,S,G,pf,gpArr3,mpArr3,weightArr3);  
            end
        end

        % Store updated policy functions  
        if strcmp(O.alg,'ART')
            pf_pigap_up(inode) = argzero(1);     
            pf_n_up(inode) = argzero(2);
            pf_q_up(inode) = argzero(3);
            pf_mc_up(inode) = argzero(4);
        elseif strcmp(O.alg,'Gust')
            pf_pigap_up(inode) = argzero(1);     
            pf_n_up(inode) = argzero(2);
            pf_n_zlb_up(inode) = argzero(3);
            pf_q_up(inode) = argzero(4);
            pf_mc_up(inode) = argzero(5);
        end
    end

    % Policy function distances
    dist_pigap = abs(pf_pigap_up - pf.pigap);
    dist_n = abs(pf_n_up - pf.n);
    dist_q = abs(pf_q_up - pf.q);
    dist_mc = abs(pf_mc_up - pf.mc);
    if strcmp(O.alg,'Gust')
        dist_n_zlb = abs(pf_n_zlb_up - pf.n_zlb);
    end

    % Maximum distance
    if strcmp(O.alg,'ART')
        dist_max = max([dist_pigap(:)',dist_n(:)',dist_q(:)',dist_mc(:)']);
    elseif strcmp(O.alg,'Gust')
        dist_max = max([dist_pigap(:)',dist_n(:)',dist_n_zlb(:)',dist_q(:)',dist_mc(:)']);
    end
    % Update policy functions
    pf.pigap = pf_pigap_up;
    pf.n = pf_n_up;
    pf.q = pf_q_up;
    pf.mc = pf_mc_up;
    if strcmp(O.alg,'Gust')
        pf.n_zlb = pf_n_zlb_up;
    end

    % Find where ZLB binds
    % Production function (2)
    y = (G.k_gr./G.g_gr).^P.alpha.*pf.n.^(1-P.alpha); 
    % Real gdp
    rgdp = G.c_gr + G.x_gr;
    rgdpp = (1-P.varphi.*(pf.pigap-1).^2/2).*y;
    % Output growth
    rgdpg = G.g_gr.*rgdpp./(P.g.*rgdp);    
    %   Interest rate rule
    inp = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi.*rgdpg.^P.phiy).^(1-P.rhoi).*exp(G.mp_gr);
    %   Percent nodes binding
    locs = find(inp <= 1);
    perbind = 100*numel(locs)/G.nodes;

    % Stopping reasons
    if dist_max > 5
        reason = 1;
    end
    
    if strcmp(O.alg,'ART')
        if (all(pf_n_up(:) < 0) || any(pf_pigap_up(:) < 0.5))
            reason = 2;
        end
    elseif strcmp(O.alg,'Gust')
        if (all(pf_n_up(:) < 0) || all(pf_n_zlb_up(:) < 0) || any(pf_pigap_up(:) < 0.5))
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
    fname = ['solution' O.it O.alg '_5'];    
    save(['solutions/' fname],'pf','O','P','S','G','V');%,'R');
end
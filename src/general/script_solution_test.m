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

%% Run Policy Function Iteration Algorithm

% Obtain Guess
pf = guess(V,P,S,G);
 
% Exogenous processes
gpArr3 = G.epsg_nodes(:,ones(O.epss_pts,1),ones(O.epsmp_pts,1));   

% Preallocate arrays to store policy function updates
%pf_c_up = zeros(G.griddim);
pf_pigap_up = zeros(G.griddim);
pf_n_up = zeros(G.griddim);
pf_q_up = zeros(G.griddim);
%pf_ups_up = zeros(G.griddim);
pf_mc_up = zeros(G.griddim);

it = 1;                                 % Iteration Counter
converged = -1;                         % Convergence Flag
reason = 0; 							% Stopping reason
dist_max = 0;                           % Max distance vector
while converged == -1
    istart = tic;                       % Iteration timer start
    for inode = 1:G.nodes
        % Find optimal policy functions on each node  
        start = [pf.pigap(inode),pf.n(inode),pf.q(inode),pf.mc(inode)]';%,pf.ups(inode)]';
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
            argzero = csolve('eqm',start,[],1e-4,10,state,...
                      O,P,S,G,pf,gpArr3,weightArr3);
        elseif strcmp(O.it,'fp')
            argzero = eqm_fp(start,state,...
                            O,P,S,G,pf,gpArr3,weightArr3);
        end
%         argzero = csolve('eqm',start,[],1e-4,10,state,...
%             O,P,S,G,pf,...
%             gpArr3,weightArr3);
        % Store updated policy functions  
        pf_pigap_up(inode) = argzero(1);     
        pf_n_up(inode) = argzero(2);
        pf_q_up(inode) = argzero(3);
        %pf_ups_up(inode) = argzero(5);
        pf_mc_up(inode) = argzero(4);
    end

    % Policy function distances
    dist_pigap = abs(pf_pigap_up - pf.pigap);
    dist_n = abs(pf_n_up - pf.n);
    dist_q = abs(pf_q_up - pf.q);
    %dist_ups = abs(pf_ups_up - pf.ups);
    dist_mc = abs(pf_mc_up - pf.mc);

    % Maximum distance
    dist_max = max([dist_pigap(:)',dist_n(:)',dist_q(:)',dist_mc(:)']);

    % Update policy functions
    pf.pigap = pf_pigap_up;
    pf.n = pf_n_up;
    pf.q = pf_q_up;
    %pf.ups = pf_ups_up;
    pf.mc = pf_mc_up;

    % Find where ZLB binds
    % Production function (2)
    y = (G.k_gr./G.g_gr).^P.alpha.*pf.n.^(1-P.alpha); 
    % Real gdp
    rgdp = G.c_gr + G.x_gr;
    rgdpp = (1-P.varphi.*(pf.pigap-1).^2/2).*y;
    % Output growth
    rgdpg = G.g_gr.*rgdpp./(P.g.*rgdp);    
    % Notional Interest Rate (9)
%    inp = in^P.rhoi*(S.i*pigap^P.phipi*rgdpg^P.phiy)^(1-P.rhoi)*exp(mp); 
%     %   HH FOC utilization (1)
%     rk = S.rk*exp(P.sigups*(pf.ups-1));
%     %   Production function (2)
%     yf = (pf.ups.*G.k_gr./G.g_gr).^P.alpha.*pf.n.^(1-P.alpha);
%     % Utilization definition (3)
%     u = S.rk*(exp(P.sigups*(pf.ups-1))-1)/P.sigups;
%     %   Firm FOC capital (4)
%     mc = rk.*pf.ups.*G.k_gr./(P.alpha*G.g_gr.*yf);
%     %   Firm FOC labor (5)
%     wp = (1-P.alpha)*mc.*yf./pf.n;
%     %   Real wage growth gap (6)
%     wg = pf.pigap.*G.g_gr.*wp./(P.g*G.w_gr);
%     %   Output definition (7)
%     yp = (1-(P.varphip*(pf.pigap-1).^2)/2-P.varphiw*(wg-1).^2/2).*yf - u.*G.k_gr./G.g_gr;
%     %   Lagged ARC (12)
%     y = G.c_gr+G.x_gr;
%     %   Output growth gap (8)
%     yg = G.g_gr.*yp./(P.g*y);
    %   Interest rate rule
    inp = G.in_gr.^P.rhoi.*(S.i*pf.pigap.^P.phipi.*rgdpg.^P.phiy).^(1-P.rhoi).*exp(G.mp_gr);
    %   Percent nodes binding
    locs = find(inp <= 1);
    perbind = 100*numel(locs)/G.nodes;

    % Stopping reasons
    if dist_max > 5
        reason = 1;
    elseif (all(pf_n_up(:) < 0) || any(pf_pigap_up(:) < 0.5))
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
    fname = ['solution' O.it];    
    save(['solutions/' fname],'pf','O','P','S','G','V');%,'R');
end
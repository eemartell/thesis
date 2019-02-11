% Solve the full nonlinear model under alternative policy parameters
clear
clc
%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------
if exist('Rules\results.mat','file')
    % Load baseline solution and structures
    load('Rules\results.mat','baselinename')
    load([baselinename '.mat'],'pf','O','P','S','G');
    % Create table of alternative parameters
    % Baseline
    baseline = [0 0.025 1.500 ];
    % Table 1: phiy
    table1 = [  1 0.000	1.500  ];  
    % Table 2: phipi                           
    table2 = [  2 0.025	1.750      
                2 0.025	2.000  ];
    tables = [baseline;table1;table2];
    % Create results matrix
    save('Rules\results.mat','tables','-append')
else
    disp('Please run script_TL_onetime.m first to get baseline solution.')
    break
end
nparams = size(tables,1);

% Solve model for each combination of phiy and phipi
perbind = 0;
itable_old = 0;
for iparam = 1:nparams
    % Set paramters
    P.phiy = tables(iparam,2);
    P.phipi = tables(iparam,3);
    suffix = ['_phiy' num2str(1e3*P.phiy) '_phipi' num2str(1e2*P.phipi)];
    savename = [baselinename suffix '.mat'];
    % Retrieve initial policy functions
    if tables(iparam,1) == 0 || tables(iparam,1) ~= itable_old
        load([baselinename '.mat'],'pf')
    end
    itable_old = tables(iparam,1);
    disp(['Solving: ' savename])
    [pf.n,pf.pi,pf.i,converged,reason,it_last,dist_last] = ...
        Fscript(O.nstates,O.npfs,O.nshocks,O.zlbflag,...
                P.alpha,P.pi,P.phipi,P.phiy,P.sigma,P.eta,P.varphi,...
                P.delta,P.nu,P.zbar,P.rhoz,P.beta,P.rhobeta,P.theta,...
                S.chi,S.r,S.y,P.tol,...
                G.k_grid,G.z_grid,G.beta_grid,...
                G.e_weight,G.e_nodes,G.u_weight,G.u_nodes, ...
                pf.n,pf.pi,pf.i);
    if converged
        save(savename,'pf','O','P','S','G')
    end
end
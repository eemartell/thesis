% Calculate simulation statistics
clear
clc
  
%--------------------------------------------------------------------------
% Run simulations
%--------------------------------------------------------------------------
results_fname = 'Rules\results';
if exist([results_fname '.mat'],'file')
    load([results_fname '.mat'])
    if exist('tables','var')
        if ~exist('sims','var')
            runsims = 'on';
        else
            runsims = 'off';
        end
    else
        disp('Run script_tables.m first to get the alternative solutions.')
        break
    end
else
    disp('Run script_TL_onetime.m first to get the baseline solution.')
    break
end

% Simulation parameters
nburn = 10000;  % Periods to get stochastic steady state
nplot = 500000;
npers = nburn+nplot;

% Set random number seed and draw shocks
mtstream = RandStream('mt19937ar','seed',2);
RandStream.setGlobalStream(mtstream);
shockse = randn(nplot,1);
shocksu = randn(nplot,1);

nparams = size(tables,1);
% Simulate for each parameterization
sims = zeros(nparams,5);
for iparam = 1:nparams
    % Load ZLB solution
    P.phiy = tables(iparam,2);
    P.phipi = tables(iparam,3);
    suffix = ['_phiy' num2str(1e3*P.phiy) '_phipi' num2str(1e2*P.phipi)];
    temp = load([baselinename suffix '.mat'],'pf','O','P','S','G');
    pf = temp.pf;
    O = temp.O;
    P = temp.P;
    S = temp.S;
    G = temp.G;
    V = variables;
    [T, M, eu] = linmodel(P,S,V);
    %--------------------------------------------------------------------------
    % Simulate
    %--------------------------------------------------------------------------
    disp(['Simulating: ' suffix]);
    % Technology shocks
    e = zeros(npers,1);
    e(nburn+1:end) = P.sige*shockse;
    % Discount factor shocks
    u = zeros(npers,1);
    u(nburn+1:end) = P.sigu*shocksu;
    % Simulate
    paths = zeros(V.nvar,npers);
    for t = 2:npers
        paths(:,t) = T*paths(:,t-1)+M*[e(t) u(t)]';
    end
    paths = paths';

    % Truncate paths
    paths = paths(nburn+1:end,:);

    % Output
    sims(iparam,:) = [P.phiy P.phipi ...
                      100*std(paths(:,V.y)) ...
                      100*std(paths(:,V.pi)) ...
                      100*std(paths(:,V.r))];
end

% Display simulation results
if exist('sims','var')
    disp('phiy          phipi    output std    inf. std      r std')
    disp(num2str(sims)) 
end
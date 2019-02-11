% Calculate simulation statistics
clear
clc
%--------------------------------------------------------------------------
% Run simulations
%--------------------------------------------------------------------------
if exist('Rules\results.mat','file')
    load('Rules\results.mat')
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

if strcmp(runsims,'on')
    nparams = size(tables,1);
    % Simulate for each parameterization
    sims = zeros(nparams,9);
    if matlabpool('size') == 0
        matlabpool
    end
    parfor iparam = 1:nparams
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
        % Production function
        y = G.z_gr.*G.k_gr.^P.alpha.*pf.n.^(1-P.alpha); 
        %   Interest rate rule     
        r = max(1,S.r*(pf.pi/P.pi).^P.phipi.*(y./O.ybench).^P.phiy);
        zlbbinds = length(find(r == 1))/G.nodes;
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
        V = variables;
        paths = simnonlin(pf,O,P,S,G,V,npers,e,u,0);

        % Truncate paths
        paths = paths(nburn+1:end,:);

        % Probability of hitting ZLB
        binds = paths(:,V.r) == 1;
        nbinds = sum(binds(:));

        % Longest ZLB event
        events = 0;
        ievent = 1; 
        temp = 0;
        for j = 2:nplot
            if binds(j) == 1 && binds(j-1) == 0
                b = j;
            elseif binds(j) == 0 && binds(j-1) == 1
                events(ievent) = j - b;
                ievent = ievent + 1;
            end
        end
        
        % Output
        sims(iparam,:) = [P.phiy P.phipi ...
                          100*zlbbinds 100*nbinds/nplot ...
                          mean(events) max(events) ...
                          100*std(paths(:,V.y))/mean(paths(:,V.y)) ...
                          100*std(paths(:,V.pi))/mean(paths(:,V.pi)) ...
                          100*std(paths(:,V.r))/mean(paths(:,V.r))];
    end
    save('Rules\results.mat','sims','-append')
else
    load('Rules\results.mat')
end

% Display simulation results
if exist('sims','var')
    disp('phiy          phipi     % of ss      % of sim     avg ZLB        max ZLB  output std    inf. std      r std')
    disp(num2str(sims))
    % Output Results to LaTeX    
    %  Comparison across phiy
    fc = {'$%.3f$','$%.2f$','$%.3f$','$%g$','$%.4f$','$%.4f$','$%.4f$'};
    matrix2latex(sims(1:6,[1,4:end]),'','formatColumns',fc)
    %  Comparison across phipi
    fc = {'$%.3f$','$%.2f$','$%.3f$','$%g$','$%.4f$','$%.4f$','$%.4f$'};
    matrix2latex(sims([1,7:11],[2,4:end]),'','formatColumns',fc)    
end
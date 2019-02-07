function sims = simulation(pf,P,S,G,V,epsg,epss,epsmp,varargin)

%----------------------------------------------------------------------
% Initialize Simulation
%----------------------------------------------------------------------
% Containers for simulated variables
[npers,nsims] = size(epsg);
sims = zeros(npers,V.nvar-V.nfore,nsims);

if isempty(varargin)
    % Initialize state variables at deterministic steady state
    sims(1,V.g,:) = P.g;
    sims(1,V.s,:) = P.s;
    sims(1,V.in,:) = S.in;
    sims(1,V.c,:) = S.c;
else
    if size(varargin{1},2) == 1
        % Initialize state variables at stochastic steady state
        sims(1,:,:) = repmat(varargin{1},[1 1 nsims]);
    elseif size(varargin{1},2) > 1
        % Initialize state variables at draw from ergodic distribution
        sims(1,:,:) = varargin{1};
    end
end

% Create 3rd dimension to match path
epsgTemp = permute(epsg,[1,3,2]);
epssTemp = permute(epss,[1,3,2]);
epsmpTemp = permute(epsmp,[1,3,2]);

%----------------------------------------------------------------------
% Simulate time series
%----------------------------------------------------------------------
for t = 2:npers
    % Growth shock
    sims(t,V.g,:) = P.g + P.sigg*epsgTemp(t,1,:);
    % Preference shock
    sims(t,V.s,:) = (1-P.rhos)*P.s + P.rhos*sims(t-1,V.s,:) + P.sigs*epssTemp(t,1,:);
    % Monetary policy shock
    sims(t,V.mp,:) = P.sigmp.*epsmpTemp(t,1,:);
    % Evaluate policy functions
    [sims(t,V.c,:),pigap(1,1,:)] = ...
        Fallterp52c_F( ...
            G.g_grid,G.s_grid,G.mp_grid,G.in_grid,G.c_grid,...
            squeeze(sims(t,V.g,:)),squeeze(sims(t,V.s,:)),squeeze(sims(t,V.mp,:)),...
                squeeze(sims(t-1,V.in,:)),squeeze(sims(t-1,V.c,:)),...
            pf.c,pf.pigap);
    sims(t,V.pi,:) = P.pi*pigap;
    % Output growth gap
    sims(t,V.yg,:) = sims(t,V.g,:).*sims(t,V.c,:)./(P.g*sims(t-1,V.c,:));
    % Interest rate rule
    sims(t,V.in,:) = sims(t-1,V.in,:).^P.rhoi.*(S.i*pigap.^P.phipi.* ...
        sims(t,V.yg,:).^P.phiy).^(1-P.rhoi).*exp(sims(t,V.mp,:));
    sims(t,V.i,:) = max(1,sims(t,V.in,:));
    % Inverse marginal utility of consumption
    sims(t,V.lam,:) = sims(t,V.c,:) - P.h*sims(t-1,V.c,:)./sims(t,V.g,:);
    % Output
    sims(t,V.y,:) = sims(t,V.c,:)./(1-P.varphip*(pigap-1).^2/2);
    % Real wage
    sims(t,V.w,:) = S.chi*sims(t,V.y,:).^P.eta.*sims(t,V.lam,:);
    
    % Check for complex
    if mod(t,10000) == 0
        if ~isreal(sims(t,:,:))
            disp('Warning: complex values in simulation');
        end
    end
end
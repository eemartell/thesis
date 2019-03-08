function sims = simulation_gustetal(pf,P,S,G,V,epsg,epss,epsmp,varargin)

%----------------------------------------------------------------------
% Initialize Simulation
%----------------------------------------------------------------------
% Containers for simulated variables
[npers,nsims] = size(epsg);
sims = zeros(npers,V.nplotvar,nsims);

if isempty(varargin)
    % Initialize state variables at deterministic steady state
    sims(1,V.g,:) = P.g;
    sims(1,V.s,:) = P.s;
    sims(1,V.in,:) = S.in;
    sims(1,V.c,:) = S.c;
    sims(1,V.k,:) = S.k;
    sims(1,V.x,:) = S.x;
    sims(1,V.y,:) = sims(1,V.c,:) + sims(1,V.x,:);
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
    % Risk Premium shock
    sims(t,V.s,:) = (1-P.rhos)*P.s + P.rhos*sims(t-1,V.s,:) + P.sigs*epssTemp(t,1,:);
    % Monetary policy shock
    sims(t,V.mp,:) = P.sigmp.*epsmpTemp(t,1,:);
    % Evaluate policy functions
    [pigap(1,1,:),sims(t,V.n,:),sims(t,V.q,:),sims(t,V.mc,:)] = ...
        Fallterp74c_F( ...
            G.g_grid,G.s_grid,G.mp_grid,...
            G.in_grid,G.c_grid,G.k_grid,G.x_grid,...
            squeeze(sims(t,V.g,:)),squeeze(sims(t,V.s,:)),...
            squeeze(sims(t,V.mp,:)),squeeze(sims(t-1,V.in,:)),...
            squeeze(sims(t-1,V.c,:)),squeeze(sims(t-1,V.k,:)),...
            squeeze(sims(t-1,V.x,:)),...
            pf.pigap,pf.n,pf.q,pf.mc);
    [pigap(1,1,:),sims(t,V.n_zlb,:),sims(t,V.q,:),sims(t,V.mc,:)] = ...
        Fallterp74c_F( ...
            G.g_grid,G.s_grid,G.mp_grid,...
            G.in_grid,G.c_grid,G.k_grid,G.x_grid,...
            squeeze(sims(t,V.g,:)),squeeze(sims(t,V.s,:)),...
            squeeze(sims(t,V.mp,:)),squeeze(sims(t-1,V.in,:)),...
            squeeze(sims(t-1,V.c,:)),squeeze(sims(t-1,V.k,:)),...
            squeeze(sims(t-1,V.x,:)),...
            pf.pigap,pf.n_zlb,pf.q,pf.mc);        
    sims(t,V.pi,:) = P.pi*pigap;
    % Production function (2)
    sims(t,V.yf,:) = (sims(t-1,V.k,:)./sims(t,V.g,:)).^P.alpha.*sims(t,V.n,:).^(1-P.alpha); 
%     % HH FOC utilization (1)
%     sims(t,V.rk,:) = P.alpha*sims(t,V.mc,:).*sims(t,V.g,:).*sims(t,V.yf,:)./sims(t-1,V.k,:); 
%     % Firm FOC labor (5)
%     sims(t,V.w,:) = (1-P.alpha)*sims(t,V.mc,:).*sims(t,V.yf,:)./sims(t,V.n,:);
%     sims(t,V.c,:) = sims(t,V.w,:)./(S.chi*sims(t,V.n,:).^P.eta) + P.h*sims(t-1,V.c,:)./sims(t,V.g,:);
    % Output definition (7)
    sims(t,V.y,:) = (1-P.varphi*(pigap-1).^2/2).*sims(t,V.yf,:);
    % Output growth gap (8)
    sims(t,V.yg,:) = sims(t,V.g,:).*sims(t,V.y,:)./(P.g*sims(t-1,V.y,:));
    % Interest rate rule (9)
    sims(t,V.in,:) = ...
        sims(t-1,V.in,:).^P.rhoi.*(S.i*pigap.^P.phipi.* ...
        sims(t,V.yg,:).^P.phiy).^(1-P.rhoi).*exp(sims(t,V.mp,:));
    sims(t,V.i,:) = max(1,sims(t,V.in,:));
    n_combined = sims(t,V.n,:).*(sims(t,V.in,:)>1) + sims(t,V.n_zlb,:).*(sims(t,V.in,:)<=1);
    % ARC (13)    c_combined = sims(t,V.c,:).*(sims(t,V.in,:)>1) + sims(t,V.c_zlb,:).*(sims(t,V.in,:)<=1); %%%combined with c's instead

%     %Recalculate (try recaluclating another time)
%     % Production function (2)
%     sims(t,V.yf,:) = (sims(t-1,V.k,:)./sims(t,V.g,:)).^P.alpha.*n_combined.^(1-P.alpha); 
%     % Output definition (7)
%     sims(t,V.y,:) = (1-P.varphi*(pigap-1).^2/2).*sims(t,V.yf,:);
%     % Output growth gap (8)
%     sims(t,V.yg,:) = sims(t,V.g,:).*sims(t,V.y,:)./(P.g*sims(t-1,V.y,:));
    
    % HH FOC utilization (1)
    sims(t,V.rk,:) = P.alpha*sims(t,V.mc,:).*sims(t,V.g,:).*sims(t,V.yf,:)./sims(t-1,V.k,:); 
    % Firm FOC labor (5)
    sims(t,V.w,:) = (1-P.alpha)*sims(t,V.mc,:).*sims(t,V.yf,:)./n_combined;
    sims(t,V.c,:) = sims(t,V.w,:)./(S.chi*n_combined.^P.eta) + P.h*sims(t-1,V.c,:)./sims(t,V.g,:);
    sims(t,V.x,:) = sims(t,V.y,:) - sims(t,V.c,:);    
    % Investment growth gap (14)
    sims(t,V.xg,:) = sims(t,V.g,:).*sims(t,V.x,:)./(P.g*sims(t-1,V.x,:));
    % Law of motion for capital (15)
    sims(t,V.k,:) = ...
        (1-P.delta)*(sims(t-1,V.k,:)./sims(t,V.g,:))+ ...
        sims(t,V.x,:).*(1-P.nu*(sims(t,V.xg,:)-1).^2/2);
    sims(t,V.lam,:) = sims(t,V.c,:) - P.h*sims(t-1,V.c,:)./sims(t,V.g,:);
%     % Check for complex
%     if mod(t,10) == 0
%         if ~isreal(sims(t,:,:))
%             disp(['Warning: complex values in simulation period ', num2str(t)]);
%         end
%     end
end
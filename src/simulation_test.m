function sims = simulation_test(pf,P,S,G,V,e,u,v,varargin)

%----------------------------------------------------------------------
% Initialize Simulation
%----------------------------------------------------------------------
[npers,nsims] = size(e);
if isempty(varargin)
    nparticles = 1;
    sims = zeros(npers,V.nplotvar,nsims);
    % States
    sims(1,V.g,:) = P.g;
    sims(1,V.s,:) = P.s;
    sims(1,V.mp,:) = 0;
    sims(1,V.in,:) = S.i;
    % Shocks
    eTemp = permute(e,[1 3 2]);
    uTemp = permute(u,[1 3 2]);
    vTemp = permute(v,[1 3 2]);
else
    statevar = varargin{1};
    [~,nparticles] = size(statevar);
    if nparticles == 1
        % States
        sims(1,:,:) = repmat(statevar,[1 nsims]);
        % Shocks
        eTemp = permute(e,[1 3 2]);
        uTemp = permute(u,[1 3 2]);
        vTemp = permute(v,[1 3 2]);
    elseif nparticles > 1
        npers = 2;
        % States
        sims(1,:,:) = statevar;
        % Shocks
        eTemp = zeros(2,1,nparticles);
        uTemp = zeros(2,1,nparticles);
        vTemp = zeros(2,1,nparticles);
        eTemp(2,1,:) = e;
        uTemp(2,1,:) = u;
        vTemp(2,1,:) = v;
    end 
end

%----------------------------------------------------------------------
% Simulate time series
%----------------------------------------------------------------------

for t = 2:npers
    % Productivity growth rate
    sims(t,V.g,:) = (1-P.rhog)*P.g + P.rhog*sims(t-1,V.g,:) + P.sige.*eTemp(t,1,:);
    % Risk premium state
    sims(t,V.s,:) = (1-P.rhos)*P.s + P.rhos*sims(t-1,V.s,:) + P.sigu.*uTemp(t,1,:);
    % Monetary policy state
    sims(t,V.mp,:) = P.rhomp*sims(t-1,V.mp,:) + P.sigv.*vTemp(t,1,:);
    % Evaluate policy functions
    [sims(t,V.c,:),pigap(1,1,:)] = Fallterp42c_F(...
        G.g_grid,G.s_grid,G.mp_grid,G.in_grid,...
        squeeze(sims(t,V.g,:)),squeeze(sims(t,V.s,:)),squeeze(sims(t,V.mp,:)),...
        squeeze(sims(t-1,V.in,:)),...
        pf.hh,pf.firm);
    sims(t,V.pi,:) = P.pi*pigap;
    % ARC
    sims(t,V.y,:) = sims(t,V.c,:)./(1-P.varphi*(pigap-1).^2/2);
    % Interest rate rule    
    sims(t,V.in,:) = sims(t-1,V.in,:).^P.rhoi.*(S.i*pigap.^P.phipi)...
        .^(1-P.rhoi).*exp(sims(t,V.mp,:));
    if P.zlbflag
        sims(t,V.i,:) = max(1,sims(t,V.in,:));
    else
        sims(t,V.i,:) = sims(t,V.in,:);
    end
    % Production function
    sims(t,V.n,:) = sims(t,V.y,:);
    % FOC Labor
    sims(t,V.w,:) = S.chi*sims(t,V.n,:).^P.eta.*(sims(t,V.c,:));
    % Lambda
    sims(t,V.lam,:) = sims(t,V.c,:); 
    % Check for complex
    if ~isreal(sims(t,:,:))
        disp('Warning: complex values in simulation');
    end
end
if nparticles > 1
    sims = squeeze(sims(2,:,:));
end

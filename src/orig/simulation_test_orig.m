function sims = simulation_test(pf,P,S,G,V,e,u,v,varargin)

%----------------------------------------------------------------------
% Initialize Simulation
%----------------------------------------------------------------------
[npers,nsims] = size(e);
if isempty(varargin) %initializes at deterministic steady state
    nparticles = 1;
    sims = zeros(npers,V.nplotvar+2,nsims);
    % States
    sims(1,V.g,:) = P.g;
    sims(1,V.a,:) = 1;
    sims(1,V.mp,:) = 0;
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
sims(1,V.c,:) = S.c;
sims(1,V.pi,:) = P.pi;
sims(1, V.i,:) = S.i;
sims(1,V.lam,:) = sims(1,V.c,:)/sims(1,V.a,:);

for t = 2:npers    
    % Productivity growth rate
    sims(t,V.g,:) = (1-P.rhog)*P.g + P.rhog*sims(t-1,V.g,:) + P.sige.*eTemp(t,1,:);
    % Preference
    sims(t,V.a,:) = (1-P.rhoa) + P.rhoa*sims(t-1,V.a,:) + P.sigu.*uTemp(t,1,:);
    % Monetary policy state
    sims(t,V.mp,:) = P.rhomp*sims(t-1,V.mp,:) + P.sigv.*vTemp(t,1,:);
    % Evaluate policy functions
    [sims(t,V.c,:),pigap(1,1,:)] = Fallterp32c_F(...
        G.g_grid,G.a_grid,G.mp_grid,...
        squeeze(sims(t,V.g,:)),squeeze(sims(t,V.a,:)),squeeze(sims(t,V.mp,:)),...
        pf.c,pf.pigap);
    sims(t,V.pi,:) = P.pi*pigap;
    % ARC
    sims(t,V.y,:) = sims(t,V.c,:)./(1-P.varphi*(pigap-1).^2/2);
    % Consumption Growth
    sims(t,V.cg,:) = sims(t,V.g,:).*sims(t,V.c,:)./sims(t-1,V.c,:); 
    % Interest rate rule    
    sims(t,V.in,:) = S.i*pigap.^P.phipi.*exp(sims(t,V.mp,:));
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
    sims(t,V.lam,:) = sims(t,V.c,:)./sims(t,V.a,:); 
    % Check for complex
    if ~isreal(sims(t,:,:))
        disp('Warning: complex values in simulation');
    end
    
    %%%%Calculating Euler equation errors 
    % SDF
    sdf = P.beta*sims(t-1,V.lam,:)./sims(t,V.lam,:);
    
    %%%Numerical integration
    %Compute all combinations of shocks
    Ebond = sdf./(sims(t,V.g,:)*pigap);
    Efp = sdf.*(pigap-1).*pigap*sims(t,V.y,:);
    
    %%%First-order conditions
    %Consumption euler equation
    sims(t,V.nplotvar+1,:) = 1 - (1+P.s)*sims(t-1,V.i,:)*Ebond/P.pi;
    sims(t,V.nplotvar+2,:) = P.varphi*(pigap-1)*pigap-(1-P.theta)- ...
        P.theta*sims(t-1,V.w,:)-P.varphi*Efp/sims(t-1,V.y,:);
    
end
if nparticles > 1
    sims = squeeze(sims(2,:,:));
end

function path = simulation(pf,P,S,G,V,e,u,varargin)

%----------------------------------------------------------------------
% Initialize Simulation
%----------------------------------------------------------------------
% Containers for simulated variables
[npers,nsims] = size(e);
path = zeros(npers,V.nvar,nsims);

if isempty(varargin)
    % Initialize state variables at deterministic steady state
    path(1,V.c,:) = S.c;
    path(1,V.k,:) = S.k;
    path(1,V.i,:) = S.i;
    path(1,V.rn,:) = S.r;
    path(1,V.g,:) = P.g;
    path(1,V.beta,:) = P.beta;
    path(1,V.rgdp,:) = S.c+S.i;
else
    % Initialize state variables at stochastic steady state
    path(1,:,:) = repmat(varargin{1},[1 1 nsims]);
end

% Create 3rd dimension to match path
eTemp = permute(e,[1 3 2]);
uTemp = permute(u,[1 3 2]);

%----------------------------------------------------------------------
% Simulate time series
%----------------------------------------------------------------------
disp('Simulating...')
for t = 2:npers
    % Technology
    path(t,V.g,:) = P.g*(path(t-1,V.g,:)/P.g).^P.rhog.*exp(eTemp(t,1,:));
    % Discount factor
    path(t,V.beta,:) = P.beta*(path(t-1,V.beta,:)/P.beta).^P.rhobeta.*exp(uTemp(t,1,:));
    % Evaluate policy functions
    [path(t,V.n,:),path(t,V.pi,:),path(t,V.psi,:)] = ...
        Fallterp63c_F(G.c_grid,G.k_grid,G.i_grid,G.rn_grid,G.g_grid,G.beta_grid,...
                    squeeze(path(t-1,V.c,:)),squeeze(path(t-1,V.k,:)),...
                    squeeze(path(t-1,V.i,:)),squeeze(path(t-1,V.rn,:)),...
                    squeeze(path(t,V.g,:)),squeeze(path(t,V.beta,:)),...
                    pf.n,pf.pi,pf.psi); 
    % Production function
    path(t,V.y,:) = (path(t-1,V.k,:)./path(t,V.g,:)).^P.alpha.*path(t,V.n,:).^(1-P.alpha); 
    % Real GDP
    path(t,V.rgdp,:) = (1-P.varphi*(path(t,V.pi,:)/P.pi-1).^2/2).*path(t,V.y,:);
    % Interest rate rule    
    path(t,V.rn,:) = path(t-1,V.rn,:).^P.rhor.*(S.r*(path(t,V.pi,:)/P.pi).^P.phipi.*...
                     (path(t,V.g,:).*path(t,V.rgdp,:)./(P.g*path(t-1,V.rgdp,:))).^P.phiy).^(1-P.rhor);
    path(t,V.r,:) = max(P.zlb,path(t,V.rn,:));
    % Firm FOC labor
    path(t,V.w,:) = (1-P.alpha)*path(t,V.psi,:).*path(t,V.y,:)./path(t,V.n,:);
    % FOC labor
    path(t,V.c,:) = path(t,V.w,:)./(S.chi*path(t,V.n,:).^P.eta)+P.h*path(t-1,V.c,:)./path(t,V.g,:);
    % Aggregate resource constraint
    path(t,V.i,:) = path(t,V.rgdp,:) - path(t,V.c,:);    
    % Investment adjustment costs
    iac = path(t,V.i,:).*path(t,V.g,:)./(path(t-1,V.i,:)*P.g);
    % Law of motion for capital
    path(t,V.k,:) = (1-P.delta)*(path(t-1,V.k,:)./path(t,V.g,:))+path(t,V.i,:).*(1-P.nu*(iac-1).^2/2);

    % Check for complex
    if mod(t,10000) == 0;
        if ~isreal(path(t,:,:))
            disp('Warning: complex values in simulation');
        end
    end
end
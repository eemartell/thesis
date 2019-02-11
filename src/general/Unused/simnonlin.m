 function path = simnonlin(pf,O,P,S,G,V,npers,e,u,expvarflag,varargin)

%----------------------------------------------------------------------
% Initialize Simulation
%----------------------------------------------------------------------
% Containers for simulated variables  
nsims = size(e,2);  
path = zeros(npers,V.nplotvar,nsims);

if isempty(varargin)
    % Initialize state variables at deterministic steady state
    path(1,V.k) = S.k;
    path(1,V.z) = P.zbar;
    path(1,V.beta,:) = P.beta;
else
    % Initialize state variables at stochastic steady state
    path(1,:,:) = repmat(varargin{1},[1 1 nsims]);
end

eTemp = permute(e,[1 3 2]);
uTemp = permute(u,[1 3 2]);
%----------------------------------------------------------------------
% Simulate time series
%----------------------------------------------------------------------
for t = 2:npers
	% Technology
    path(t,V.z,:) = P.zbar*(path(t-1,V.z,:)/P.zbar).^P.rhoz.*exp(eTemp(t,1,:));
    % Discount factor
    path(t,V.beta,:) = P.beta*(path(t-1,V.beta,:)/P.beta).^P.rhobeta.*exp(uTemp(t,1,:));
    % Evaluate policy functions
    [path(t,V.n,:),path(t,V.pi,:),path(t,V.i,:)] = ...
        Fallterp33c_F(G.k_grid,G.z_grid,G.beta_grid, ...
                      path(t-1,V.k,:),path(t,V.z,:),path(t,V.beta,:), ...
                      pf.n,pf.pi,pf.i);       
    % Production function
    path(t,V.y,:) = path(t,V.z,:).*path(t-1,V.k,:).^P.alpha.*path(t,V.n,:).^(1-P.alpha);    
    % Interest rate rule    
    path(t,V.r,:) = S.r*(path(t,V.pi,:)/P.pi).^P.phipi.*(path(t,V.y,:)./S.y).^P.phiy;
    if O.zlbflag == 1
        path(t,V.r,:) = max(1,path(t,V.r,:));
    end
    % Tobin's q
    path(t,V.q,:) = 1+P.nu*(path(t,V.i,:)./path(t-1,V.k,:)-P.delta);
    % Aggregate resource constraint
    ytil = path(t,V.y,:).*(1-(P.varphi*(path(t,V.pi,:)/P.pi-1).^2)/2);
    kac = P.nu*(path(t,V.i,:)./path(t-1,V.k,:)-P.delta).^2/2;
    path(t,V.c,:) = ytil - path(t,V.i,:) - kac.*path(t-1,V.k,:);
    % FOC Labor
    path(t,V.w,:) = S.chi*path(t,V.n,:).^P.eta.*path(t,V.c,:).^P.sigma;
    % Firm FOC labor
    path(t,V.psi,:) = path(t,V.w,:).*path(t,V.n,:)./((1-P.alpha)*path(t,V.y,:));
    % Firm FOC capital
    path(t,V.rk,:) = P.alpha*path(t,V.psi,:).*path(t,V.y,:)./path(t-1,V.k,:);    
    if ~isreal(path(t,:,:))
        display('warning: complex path -> excessive extrapolation of pfs')
    end
    % Expected variables
    if expvarflag == 1
        k = squeeze(path(t-1,V.k,:));
        z = squeeze(path(t,V.z,:));
        beta = squeeze(path(t,V.beta,:));
        i = path(t,V.i,:);
        [Einvpi,path(t,V.Ec,:),path(t,V.Er,:),path(t,V.Erk,:)] = ...
                        expvar(k,z,beta,i,...
                               P.alpha,P.pi,P.eta,P.varphi,P.phipi,P.phiy,...  
                               P.delta,P.nu,P.zbar,P.rhoz,P.beta,P.rhobeta,...
                               S.chi,S.r,S.y,...       
                               G.e_weight,G.e_nodes,G.u_weight,G.u_nodes,...
                               G.k_grid,G.z_grid,G.beta_grid,...
                               pf.n,pf.pi,pf.i);
       path(t,V.realr,:) = squeeze(path(t,V.r,:)).*Einvpi;
    end
    % Investment
    path(t,V.k,:) = path(t,V.i,:) + (1-P.delta)*path(t-1,V.k,:); 
end
% Rotemberg adjustment cost
path(:,V.rotcost,:) = (P.varphi*(path(:,V.pi,:)/P.pi-1).^2)/2;
% Adjusted output
path(:,V.ytil,:) = path(:,V.y,:).*(1-(P.varphi*(path(:,V.pi,:)/P.pi-1).^2)/2);
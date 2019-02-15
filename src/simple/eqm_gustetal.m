function R = eqm_gustetal(x,state,O,P,S,G,pf,gpArr3,mpArr3,weightArr3,varargin)
% Get original grids and GH nodes
if ~isempty(varargin)
    EEflag = 1;
    GH = varargin{1};
else
    EEflag = 0;
end

% Preallocate function output
R = zeros(size(x));
Rdim = size(R,2);

% State Values
g = state(1);       %Growth state current period
s = state(2);       %Preference state current period
mp = state(3);      %Monetary policy state current period
in = state(4);  	%Notional interest rate last period

for icol = 1:Rdim   
% Policy Function Guesses
cp = x(1,icol);      %Consumption policy current period, non-ZLB
%cp_zlb = cp;
cp_zlb = x(2,icol);  %Consumption policy current period, ZLB
pigap = x(3,icol);   %Inflation gap policy current period    
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
% Interest rate rule (2,3)
inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
i = max(1,inp);
if inp > 1
    lam = cp;
    % Aggregate resource constraint (6)    
    y = cp/(1-P.varphi*(pigap-1)^2/2);
else
    lam = cp_zlb;
    % Aggregate resource constraint (6)    
    y = cp_zlb/(1-P.varphi*(pigap-1)^2/2); 
end
% FOC Labor (7)
w = S.chi*y^P.eta*lam;
%----------------------------------------------------------------------
% Linear interpolation of the policy variables
%----------------------------------------------------------------------
if ~EEflag         
    [cppArr3,pigappArr3] = Fallterp423_R(...
        O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
        G.in_grid,...
        inp,...
        pf.c,pf.pigap);
    [cppArr3_zlb,pigappArr3] = Fallterp423_R(...
        O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
        G.in_grid,...
        inp,...
        pf.c_zlb,pf.pigap);
else
    % Growth rate
    gpVec = (1-P.rhog)*P.g + P.rhog*g + GH.e_nodes;
    % Risk premium 
    spVec = (1-P.rhos)*P.s + P.rhos*s + GH.u_nodes;
    % Monetary policy
    mpVec = GH.v_nodes;    
    [cppArr3,pigappArr3] = allterp423(...
                            G.g_grid,G.s_grid,G.mp_grid,G.in_grid,...
                            gpVec,spVec,GH.v_nodes,inp,...
                            pf.c,pf.pigap);
    [cppArr3_zlb,pigappArr3] = allterp423(...
                            G.g_grid,G.s_grid,G.mp_grid,G.in_grid,...
                            gpVec,spVec,GH.v_nodes,inp,...
                            pf.c_zlb,pf.pigap);                  
    gpArr3 = gpVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
    mpArr3 = mpVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
end

%----------------------------------------------------------------------        
% Solve for variables inside expectations
%---------------------------------------------------------------------- 
% Back out pigap 
inpArr3 = inp^P.rhoi.*(S.i*pigappArr3.^P.phipi).^(1-P.rhoi).*exp(mpArr3); %(2)
cppArr3_combined = cppArr3.*(inpArr3>1) + cppArr3_zlb.*(inpArr3<=1);
% Aggregate resource constraint (6)  
ypArr3 = cppArr3_combined./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor
lampArr3 = cppArr3_combined;
sdfArr3 = P.beta*lam./lampArr3;
%----------------------------------------------------------------------
% Numerical integration
%----------------------------------------------------------------------
% Compute all combinations of shocks
EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*pigappArr3);
EfpArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*ypArr3;
% Integrate
Ebond = sum(EbondArr3(:));
Efp = sum(EfpArr3(:));
%----------------------------------------------------------------------
% First-order conditions
%----------------------------------------------------------------------
% Consumption Euler equation (3)
R(1,icol) = 1 - s*i*Ebond/P.pi;
R(2,icol) = 1 - s*Ebond/P.pi;
% Firm Pricing (Philips Curve, 5)
R(3,icol) = P.varphi*(pigap-1)*pigap-(1-P.theta)-P.theta*w-P.varphi*Efp/y;
end
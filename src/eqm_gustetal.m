function R = eqm_gustetal(x,state,O,P,S,G,pf,gpArr3,mpArr3,weightArr3,varargin)
%%%change in the stye of eqm_fp_gustetal, noting that you need residuals
%%%rather than solving for c and pi
% Get original grids and GH nodes
if ~isempty(varargin)
    EEflag = 1;
    GH = varargin{1};
else
    EEflag = 0;
end

% Preallocate function output
R = zeros(size(x));

% State Values
g = state(1);       %Growth state current period
s = state(2);       %Preference state current period
mp = state(3);      %Monetary policy state current period
in = state(4);  	%Notional interest rate last period

% Policy Function Guesses
% Put V functions here instead
Vlambdap = x(1);      %Vlambda policy current period, non-ZLB
Vlambdap_zlb = x(2);  %Vlambda policy current period, ZLB
Vpip = x(3);          % Vpi policy current period    
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
% Back out pigap
pigap = (1+sqrt((P.varphi + 4*Vpip)/P.varphi))/2; %(1)
% Interest rate rule (2,3)
inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
i = max(1,inp);
if inp > 1
    lam = 1/Vlambdap; %(4)
else
    lam = 1/Vlambdap_zlb; %(4)
end
c = lam; %(5)
% Aggregate resource constraint (6)    
y = c/(1-P.varphi*(pigap-1)^2/2);
% FOC Labor (7)
w = S.chi*y^P.eta*lam;
%----------------------------------------------------------------------
% Linear interpolation of the policy variables
%----------------------------------------------------------------------
if ~EEflag         
    [VlambdapArr3,VpipArr3] = Fallterp423_R(...
        O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
        G.in_grid,...
        inp,...
        pf.hh,pf.firm);
    [VlambdapArr3_zlb,VpipArr3] = Fallterp423_R(...
        O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
        G.in_grid,...
        inp,...
        pf.hh_zlb,pf.firm);
else
    % Growth rate
    gpVec = (1-P.rhog)*P.g + P.rhog*g + GH.e_nodes;
    % Risk premium 
    spVec = (1-P.rhos)*P.s + P.rhos*s + GH.u_nodes;
    % Monetary policy
    mpVec = GH.v_nodes;    
    [VlambdapArr3,VpipArr3] = allterp423(...
                            G.g_grid,G.s_grid,G.mp_grid,G.in_grid,...
                            gpVec,spVec,GH.v_nodes,inp,...
                            pf.hh,pf.firm);
    [VlambdapArr3_zlb,VpipArr3] = allterp423(...
                            G.g_grid,G.s_grid,G.mp_grid,G.in_grid,...
                            gpVec,spVec,GH.v_nodes,inp,...
                            pf.hh_zlb,pf.firm);                      
    gpArr3 = gpVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
    mpArr3 = mpVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
end

%----------------------------------------------------------------------        
% Solve for variables inside expectations
%---------------------------------------------------------------------- 
% Back out pigap 
pigappArr3 = (1+sqrt((P.varphi + 4*VpipArr3)/P.varphi))/2; %(1)
inpArr3 = inp^P.rhoi.*(S.i*pigappArr3.^P.phipi).^(1-P.rhoi).*exp(mpArr3); %(2)
VlambdapArr3_combined = VlambdapArr3.*(inpArr3>1) + VlambdapArr3_zlb.*(inpArr3<=1);
% Solve for variables using combined Vlambda
lampArr3 = 1/VlambdapArr3_combined; %(4)
cppArr3 = lampArr3; %(5)
% Aggregate resource constraint (6)  
ypArr3 = cppArr3./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor
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
% probably doesn't change unless they are writtin in terms of V somehow
% Consumption Euler equation (3)
R(1) = 1 - s*i*Ebond/P.pi;
R(2) = 1 - s*Ebond/P.pi;
% Firm Pricing (Philips Curve, 5)
R(3) = P.varphi*(pigap-1)*pigap-(1-P.theta)-P.theta*w-P.varphi*Efp/y;

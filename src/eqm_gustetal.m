function R = eqm_gustetal(x,state,O,P,S,G,pf,gpArr3,apArr3,mpArr3,weightArr3,varargin)

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
g = state(1);      %Notional interest rate last period
a = state(2);       %Growth state current period
mp = state(3);       %Preference state current period
in = state(4);      %Monetary policy state current period

% Policy Function Guesses
cp = x(1);     %Consumption policy current period 
cp_zlb = x(2);     %Consumption policy current period 
pigap = x(3);	%Inflation policy current period     
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
i = max(1,inp);
if inp>1
    % Aggregate resource constraint (6)    
    y = cp/(1-P.varphi*(pigap-1)^2/2);
    % Interest rate rule (8)
    % FOC Labor (1,2)
    lam = cp/a;
else
% Aggregate resource constraint (6)    
    y = cp_zlb/(1-P.varphi*(pigap-1)^2/2);
    % Interest rate rule (8)
    % FOC Labor (1,2)
    lam = cp_zlb/a;
end
w = S.chi*a*y^P.eta*lam;
%----------------------------------------------------------------------
% Linear interpolation of the policy variables
%----------------------------------------------------------------------
%cppArr3 = pf.c;
%pigappArr3 = pf.pigap;
[cppArr3,pigappArr3] = Fallterp423_R(... %Fallterp423r???
    O.g_pts,O.a_pts,O.mp_pts,O.in_pts,...
    G.in_grid,...
    inp,...
    pf.c,pf.pigap);
[cppArr3_zlb,pigappArr3] = Fallterp423_R(... %Fallterp423r???
    O.g_pts,O.a_pts,O.mp_pts,O.in_pts,...
    G.in_grid,...
    inp,...
    pf.c_zlb,pf.pigap);
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
inpArr3 = inp^P.rhoi*(S.i*pigappArr3.^P.phipi).^(1-P.rhoi).*exp(mpArr3);
cppArr3_combined = cppArr3.*(inpArr3>1) + cppArr3_zlb.*(inpArr3<=1);
% Aggregate resource constraint  
ypArr3 = cppArr3_combined./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor
lampArr3 = cppArr3_combined./apArr3;
sdfArr3 = P.beta*lam./lampArr3;
%----------------------------------------------------------------------
% Numerical integration
%----------------------------------------------------------------------
% Weight realizations
EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*pigappArr3);
EfpArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*ypArr3;
% Integrate
Ebond = sum(EbondArr3(:));
Efp = sum(EfpArr3(:));
%----------------------------------------------------------------------
% First-order conditions
%----------------------------------------------------------------------
% Consumption Euler equation (3)
R(1) = 1 - (1+P.s)*i*Ebond/P.pi;
R(2) = 1 - (1+P.s)*Ebond/P.pi;
% Firm Pricing (Philips Curve, 5)
R(3) = P.varphi*(pigap-1)*pigap-(1-P.theta)-P.theta*w-P.varphi*Efp/y;

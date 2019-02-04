function x_up = eqm_fp(x,state,O,P,S,G,pf,gpArr3,apArr3,weightArr3)

% State Values
g = state(1);      %Notional interest rate last period
a = state(2);       %Growth state current period
mp = state(3);       %Preference state current period
in = state(4);      %Monetary policy state current period

% Policy Function Guesses
cp = x(1);     %Consumption policy current period 
pigap = x(2);	%Inflation policy current period     
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
% Aggregate resource constraint (6)    
y = cp/(1-P.varphi*(pigap-1)^2/2);
% Interest rate rule (8)
inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
if P.zlbflag
    i = max(1,inp);
else
    i = inp;
end
% FOC Labor (1,2)
lam = cp/a;
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
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
% Aggregate resource constraint  
ypArr3 = cppArr3./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor
lampArr3 = cppArr3./apArr3;
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
% Integrate
RHS_firm = 1 - P.theta + P.theta*w + P.varphi*Efp/y;
%----------------------------------------------------------------------
% First-order conditions
%----------------------------------------------------------------------
% Consumption Euler Equation
x_up(1) = P.pi*lam*a/((1+P.s)*i*Ebond);
% Firm Pricing 
x_up(2) = (1+sqrt((P.varphi+4*RHS_firm)/P.varphi))/2;
end

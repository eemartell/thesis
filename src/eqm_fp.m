function x_up = eqm_fp(x,state,O,P,S,G,pf,gpArr3,spArr3,weightArr3)

% State Values
g = state(1);      %Notional interest rate last period
s = state(2);       %Growth state current period
mp = state(3);       %Preference state current period
in = state(4);      %Monetary policy state current period

% Policy Function Guesses
cp = x(1);     %Consumption policy current period 
pigap = x(2);	%Inflation policy current period     
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
% Aggregate resource constraint (1)    
y = cp/(1-P.varphi*(pigap-1)^2/2);
% Interest rate rule (2,3)
inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
if P.zlbflag
    i = max(1,inp);
else
    i = inp;
end
% FOC Labor (4,5)
lam = cp;
w = S.chi*y^P.eta*lam;
%----------------------------------------------------------------------
% Linear interpolation of the policy variables
%----------------------------------------------------------------------
[cppArr3,pigappArr3] = Fallterp423_R(...
    O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
    G.in_grid,...
    inp,...
    pf.hh,pf.firm); %pf.c, pf.pigap
    %pf.c,pf.pigap);
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
% Aggregate resource constraint (1)  
ypArr3 = cppArr3./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor (4,5)
lampArr3 = cppArr3;
sdfArr3 = P.beta*lam./lampArr3;
%----------------------------------------------------------------------
% Numerical integration
%----------------------------------------------------------------------
% Weight realizations
EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*pigappArr3);
EfpArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*ypArr3;
% Integrate
Vlambda = s*i*sum(EbondArr3(:))/(P.pi*lam);  %%%just have as x_up(1)
Vpi = 1 - P.theta + P.theta*w + P.varphi*sum(EfpArr3(:))/y; %%%just have as x_up(2)
%----------------------------------------------------------------------
% First-order conditions
%----------------------------------------------------------------------
% Consumption Euler Equation (6)
x_up(1) = 1/Vlambda; %%%no Vlambda
% Firm Pricing (7)
x_up(2) = (1+sqrt((P.varphi+4*Vpi)/P.varphi))/2; %%%no Vpi
end

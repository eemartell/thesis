function x_up = eqm_fp_gustetal(x,state,O,P,S,G,pf,gpArr3,spArr3,weightArr3)

% State Values
g = state(1);      %Growth state current period
s = state(2);      %Risk premium state current period
mp = state(3);     %Monetary policy state current period  
in = state(4);     %Notional interest rate last period
 
%%% Unpack all 4 policy functions
% Policy Function Guesses
Vlambdap = x(1); %Vlambda policy current period
Vpip = x(2);     % Vpi policy current period   
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
%%% Solve for interest rate using non-ZLB Vs
% Back out pigap
lam = 1/Vlambdap; %(1)
c = lam; %(2)
pigap = (1+sqrt((P.varphi + 4*Vpip)/P.varphi))/2; %(3)
% Aggregate resource constraint (4)    
y = c/(1-P.varphi*(pigap-1)^2/2);
% Interest rate rule (5,6)
inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
if P.zlbflag
    i = max(1,inp);
else
    i = inp;
end
%%% If interest rate is greater than 1, solve for the rest of time t variable
%%%s using non-ZLB Vs. If not, solve for all of the time t variables using
%%%ZLB Vs.
% FOC Labor (7)
w = S.chi*y^P.eta*lam;
%----------------------------------------------------------------------
% Linear interpolation of the policy variables
%----------------------------------------------------------------------
[VlambdapArr3,VpipArr3] = Fallterp423_R(... %%%Interpolate all 4 (will need new interpolation fxn)
    O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
    G.in_grid,...
    inp,...
    pf.hh,pf.firm); %%%Interpolate all 4
    %pf.c,pf.pigap);
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
%%%solve for updated interest rate to see at which nodes ZLB binds. Construct
%%%new Vs that are non-ZLB where the interest rate is greater than 1 and ZLB
%%%otherwise. 
% Back out pigap 
lampArr3 = 1/VlambdapArr3; %(1)
cppArr3 = lampArr3; %(2)
pigappArr3 = (1+sqrt((P.varphi + 4*VpipArr3)/P.varphi))/2; %(3)
%%% Solve for all time t variables (inc ones above) using new V. 
% Aggregate resource constraint (4)  
ypArr3 = cppArr3./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor
sdfArr3 = P.beta*lam./lampArr3;
%----------------------------------------------------------------------
% Numerical integration
%----------------------------------------------------------------------
%%% Integrate with new Vs
% Weight realizations
EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*pigappArr3);
EfpArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*ypArr3;
% Integrate
x_up(1) = s*i*sum(EbondArr3(:))/(P.pi*lam); %(8)
x_up(2) = 1 - P.theta + P.theta*w + P.varphi*sum(EfpArr3(:))/y; %(9)
end

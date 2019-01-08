function x_up = eqm_fp_gustetal(x,state,O,P,S,G,pf,gpArr3,mpArr3,weightArr3)

% State Values
g = state(1);      %Growth state current period
s = state(2);      %Risk premium state current period
mp = state(3);     %Monetary policy state current period  
in = state(4);     %Notional interest rate last period
 
%%% Unpack all 4 policy functions
% Policy Function Guesses
Vlambdap = x(1);     %Vlambda policy current period
Vlambdap_zlb = x(2); 
Vpip = x(3);         % Vpi policy current period
Vpip_zlb = x(4);
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
%%% If interest rate is greater than 1, solve for the rest of time t variable
%%%s using non-ZLB Vs. If not, solve for all of the time t variables using
%%%ZLB Vs.
%%% Solve for interest rate using non-ZLB Vs
% Back out pigap
pigap = (1+sqrt((P.varphi + 4*Vpip)/P.varphi))/2; %(3)
% Interest rate rule (5,6)
inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
i = max(1,inp);
if inp > 1
    lam = 1/Vlambdap; %(1)
else
    lam = 1/Vlambdap_zlb; %(1)
    pigap = (1+sqrt((P.varphi + 4*Vpip_zlb)/P.varphi))/2; %(3)
    % Interest rate rule (5,6)
    inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
    i = max(1,inp);
end
c = lam; %(2)
% Aggregate resource constraint (4)    
y = c/(1-P.varphi*(pigap-1)^2/2);
% FOC Labor (7)
w = S.chi*y^P.eta*lam;
%----------------------------------------------------------------------
% Linear interpolation of the policy variables
%----------------------------------------------------------------------
% [VlambdapArr3,VlambdapArr3_zlb,VpipArr3,VpipArr3_zlb] = Fallterp443_R(... %%%Interpolate all 4 (will need new interpolation fxn)
%     O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
%     G.in_grid,...
%     inp,...
%     pf.hh,pf.hh_zlb,pf.firm,pf.firm_zlb); %%%Interpolate all 4
[VlambdapArr3,VpipArr3] = Fallterp423_R(...
    O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
    G.in_grid,...
    inp,...
    pf.hh,pf.firm);
[VlambdapArr3_zlb,VpipArr3_zlb] = Fallterp423_R(...
    O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
    G.in_grid,...
    inp,...
    pf.hh_zlb,pf.firm_zlb);
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
%%%solve for updated interest rate to see at which nodes ZLB binds. Construct
%%%new Vs that are non-ZLB where the interest rate is greater than 1 and ZLB
%%%otherwise. 
% Back out pigap 
pigappArr3 = (1+sqrt((P.varphi + 4*VpipArr3)/P.varphi))/2; %(3)
inpArr3 = inp^P.rhoi.*(S.i*pigappArr3.^P.phipi).^(1-P.rhoi).*exp(mpArr3);
VlambdapArr3_combined = VlambdapArr3.*(inpArr3>1) + VlambdapArr3_zlb.*(inpArr3<=1);
VpipArr3_combined = VpipArr3.*(inpArr3>1) + VpipArr3_zlb.*(inpArr3<=1);
%%% Solve for all time t variables (inc ones above) using new V.
lampArr3 = 1/VlambdapArr3_combined; %(1)
cppArr3 = lampArr3; %(2)
pigappArr3 = (1+sqrt((P.varphi + 4*VpipArr3_combined)/P.varphi))/2; %(3)
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
x_up(2) = s*sum(EbondArr3(:))/(P.pi*lam);
x_up(3) = 1 - P.theta + P.theta*w + P.varphi*sum(EfpArr3(:))/y; %(9)
x_up(4) = x_up(3); %Don't separate for Vpip_zlb???
end

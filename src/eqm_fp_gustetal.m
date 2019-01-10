function x_up = eqm_fp_gustetal(x,state,O,P,S,G,pf,gpArr3,mpArr3,weightArr3)

% State Values
g = state(1);      %Growth state current period
s = state(2);      %Risk premium state current period
mp = state(3);     %Monetary policy state current period  
in = state(4);     %Notional interest rate last period
 
% Policy Function Guesses
%%%these will be c, c_zlb, and pigap
Vlambdap = x(1);     %non-ZLB Vlambda policy current period
Vlambdap_zlb = x(2); %ZLB Vlambda policy current period
Vpip = x(3);         % Vpi policy current period
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
%%%unpacking and solving for variables will be different here
%%%solve immediately for interest rate
%%%only solve for what's needed in last steps
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
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
% Back out pigap 
%%%can solve immediately for interest rate and then split up c into regimes
pigappArr3 = (1+sqrt((P.varphi + 4*VpipArr3)/P.varphi))/2; %(1)
inpArr3 = inp^P.rhoi.*(S.i*pigappArr3.^P.phipi).^(1-P.rhoi).*exp(mpArr3); %(2)
VlambdapArr3_combined = VlambdapArr3.*(inpArr3>1) + VlambdapArr3_zlb.*(inpArr3<=1);
% Solve for variables using combined Vlambda
%%%solve for variables using combined c_combined
lampArr3 = 1/VlambdapArr3_combined; %(4)
cppArr3 = lampArr3; %(5)
% Aggregate resource constraint (6)  
ypArr3 = cppArr3./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor
sdfArr3 = P.beta*lam./lampArr3;
%----------------------------------------------------------------------
% Numerical integration
%----------------------------------------------------------------------
% Weight realizations
EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*pigappArr3);
EfpArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*ypArr3;
% Integrate
%%%solve for c and pigap here
%%%start fixing from bottom and work up so you know what you'll need
x_up(1) = s*i*sum(EbondArr3(:))/(P.pi*lam); %(8)
x_up(2) = s*sum(EbondArr3(:))/(P.pi*lam);
x_up(3) = 1 - P.theta + P.theta*w + P.varphi*sum(EfpArr3(:))/y; %(9)
end

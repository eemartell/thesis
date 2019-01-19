function x_up = eqm_fp_gustetal(x,state,O,P,S,G,pf,gpArr3,mpArr3,weightArr3)

% State Values
g = state(1);      %Growth state current period
s = state(2);      %Risk premium state current period
mp = state(3);     %Monetary policy state current period  
in = state(4);     %Notional interest rate last period
 
% Policy Function Guesses
%%%these will be c, c_zlb, and pigap
cp = x(1);      %Non-ZLB consumption policy current period 
cp_zlb = x(2);  %ZLB consumption policy current period 
pigap = x(3);   % Inflation gap policy current period
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
%%%unpacking and solving for variables will be different here
%%%solve immediately for interest rate
%%%only solve for what's needed in last steps
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
[cppArr3,pigappArr3] = Fallterp423_R(...
    O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
    G.in_grid,...
    inp,...
    pf.c,pf.pigap); %pf.c, pf.pigap
[cppArr3_zlb,pigappArr3] = Fallterp423_R(...
    O.g_pts,O.s_pts,O.mp_pts,O.in_pts,...
    G.in_grid,...
    inp,...
    pf.c_zlb,pf.pigap); %pf.c, pf.pigap
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
% Back out pigap 
%%%can solve immediately for interest rate and then split up c into regimes
inpArr3 = inp^P.rhoi.*(S.i*pigappArr3.^P.phipi).^(1-P.rhoi).*exp(mpArr3); %(2)
cppArr3_combined = cppArr3.*(inpArr3>1) + cppArr3_zlb.*(inpArr3<=1);
% Solve for variables using combined Vlambda
%%%solve for variables using combined c_combined
% Aggregate resource constraint (6)  
ypArr3 = cppArr3_combined./(1-P.varphi*(pigappArr3-1).^2/2);
% Stochastic discount factor
lampArr3 = cppArr3_combined;
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
Vlam = s*i*sum(EbondArr3(:))/(P.pi*lam); %(8) %%%just do x_up(1)
Vlam_zlb = s*sum(EbondArr3(:))/(P.pi*lam); %%%just do x_up(2)
Vpi = 1 - P.theta + P.theta*w + P.varphi*sum(EfpArr3(:))/y; %(9) %%%just do x_up(3)
%----------------------------------------------------------------------
% First-order conditions
%----------------------------------------------------------------------
% Consumption Euler Equation (6)
x_up(1) = 1/Vlam; %%%remove
x_up(2) = 1/Vlam_zlb;
% Firm Pricing (7)
x_up(3) = (1+sqrt((P.varphi+4*Vpi)/P.varphi))/2;
end

function x_up = eqm_fp_gustetal(pf0,state,O,P,S,G,pf,gpArr3,mpArr3,weightArr3)

% State Values
g = state(1);           %Growth state current period
s = state(2);           %Risk premium state current period
mp = state(3);          %Monetary policy shock current period
in = state(4);          %Notional interest rate last period
c = state(5);           %Consumption last period
k = state(6);           %Capital last period
x = state(7);           %Investment last period   

% Policy Function Guesses
pigap = pf0(1);   %Inflation gap current period  
n = pf0(2);       %Labor current period 
n_zlb = pf0(3);
q = pf0(4);       %Tobin's q current period
mc = pf0(5);      %Marginal cost current period    
%----------------------------------------------------------------------
% Current period
%----------------------------------------------------------------------
% Production function (2)
yf = (k/g)^P.alpha*n^(1-P.alpha);
% Real gdp
y = c + x;
yp = (1-P.varphi*(pigap-1)^2/2)*yf;
% Output growth
yg = g*yp/(P.g*y);    
% Notional Interest Rate (9)
inp = in^P.rhoi*(S.i*pigap^P.phipi*yg^P.phiy)^(1-P.rhoi)*exp(mp); 
% Nominal Interest Rate (10)
i = inp; 
if inp > 1
    % Firm FOC labor (5)
    w = (1-P.alpha)*mc*yf/n;
    % FOC labor
    cp = w/(S.chi*n^P.eta)+P.h*c/g;
else
    % Firm FOC labor (5)
    w = (1-P.alpha)*mc*yf/n_zlb;
    % FOC labor
    cp = w/(S.chi*n_zlb^P.eta)+P.h*c/g;
end
% Aggregate resource constraint
xp = yp - cp;
% Investment growth gap (14)
xg = g*xp/(P.g*x);    
% Law of motion for capital (15)
kp = (1-P.delta)*(k/g)+xp*(1-P.nu*(xg-1)^2/2);       
% Inverse MUC (11)
lam = cp-P.h*c/g;
%----------------------------------------------------------------------
% Linear interpolation of the policy functions 
%---------------------------------------------------------------------- 
[pigappArr3,npArr3,qpArr3,mcpArr3] = Fallterp743_R(...
    O.g_pts,O.s_pts,O.mp_pts,...
    O.in_pts,O.c_pts,O.k_pts,O.x_pts,...
    G.in_grid,G.c_grid,G.k_grid,G.x_grid,...
    inp,cp,kp,xp,...
    pf.pigap,pf.n,pf.q,pf.mc);
[~,npArr3_zlb,~,~] = Fallterp743_R(...
    O.g_pts,O.s_pts,O.mp_pts,...
    O.in_pts,O.c_pts,O.k_pts,O.x_pts,...
    G.in_grid,G.c_grid,G.k_grid,G.x_grid,...
    inp,cp,kp,xp,...
    pf.pigap,pf.n_zlb,pf.q,pf.mc);
%----------------------------------------------------------------------        
% Next period
%----------------------------------------------------------------------  
% Production function (2)
yfpArr3 = (kp./gpArr3).^P.alpha.*npArr3.^(1-P.alpha);
% Real gdp
yp = cp + xp;
yppArr3 = (1-P.varphi.*(pigappArr3-1).^2/2).*yfpArr3;
% Output growth
ygpArr3 = gpArr3.*yppArr3./(P.g.*yp);
% Notional Interest Rate (9)
inpArr3 = inp.^P.rhoi.*(S.i.*pigappArr3.^P.phipi.*ygpArr3.^P.phiy).^(1-P.rhoi).*exp(mpArr3);
npArr3_combined = npArr3.*(inpArr3>1) + npArr3_zlb.*(inpArr3<=1);
% Firm FOC capital (4)
rkpArr3 = P.alpha.*mcpArr3.*gpArr3.*yfpArr3/kp;
% Firm FOC labor (5)
wpArr3 = (1-P.alpha)*mcpArr3.*yfpArr3./npArr3_combined;
% FOC labor
cppArr3 = wpArr3./(S.chi*npArr3_combined.^P.eta)+P.h*cp./gpArr3;
% ARC
xppArr3 = yppArr3-cppArr3;
% Inverse MUC
lampArr3 = cppArr3-P.h*cp./gpArr3;
% Investment growth gap (14)
xgpArr3 = gpArr3.*xppArr3/(P.g*xp);
% Stochastic discount factor
sdfArr3 = P.beta*lam./lampArr3;
%----------------------------------------------------------------------
% Expectations
%----------------------------------------------------------------------
EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*(P.pi*pigappArr3));
EcapArr3 = weightArr3.*sdfArr3.*(rkpArr3+(1-P.delta)*qpArr3)./gpArr3;
EinvArr3 = weightArr3.*sdfArr3.*qpArr3.*xgpArr3.^2.*(xgpArr3-1)./gpArr3;
EppcArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*(yfpArr3/yf);
Ebond = sum(EbondArr3(:));
Ecap = sum(EcapArr3(:));
Einv = sum(EinvArr3(:));
Eppc = sum(EppcArr3(:));
%----------------------------------------------------------------------
% Fixed Point
%----------------------------------------------------------------------
x_up(4) = Ecap; %pf.q
pf_lam = 1/(s*i*Ebond/lam);
c_pf = pf_lam + P.h*c/g;
var = 1-(1-P.nu*P.g*Einv)/x_up(4);
xg_pf = 1/3*(2+sqrt(P.nu*(P.nu+6*var))/P.nu);
x_pf = xg_pf*P.g*x/g;
ygdp_pf = c_pf + x_pf;
RHS_firm = 1 - P.theta + P.theta*mc + P.varphi*Eppc;
x_up(1) = (1+sqrt((P.varphi+4*RHS_firm)/P.varphi))/2;

y_pf = ygdp_pf/(1-P.varphi*(x_up(1)-1)^2/2);
x_up(2) = (y_pf/(k/g)^P.alpha)^(1/(1-P.alpha)); %pf.n

x_up(5) = (w*x_up(2))/((1-P.alpha)*y_pf);

pf_lam_zlb = 1/(s*Ebond/lam);
c_pf = pf_lam_zlb + P.h*c/g;
ygdp_pf = c_pf + x_pf;

y_pf = ygdp_pf/(1-P.varphi*(x_up(1)-1)^2/2);
x_up(3) = (y_pf/(k/g)^P.alpha)^(1/(1-P.alpha)); %pf.n

end

function Res = eqm(pf0,state,O,P,S,G,pf,gpArr3,weightArr3)

% Preallocate function output
Res = zeros(size(pf0));
ncol = size(Res,2);

% State Values
g = state(1);           %Growth state current period
s = state(2);           %Risk premium state current period
mp = state(3);          %Monetary policy shock current period
in = state(4);          %Notional interest rate last period
c = state(5);           %Consumption last period
k = state(6);           %Capital last period
x = state(7);           %Investment last period   
%w = state(8);           %Real wage last period

for icol = 1:ncol
    % Policy Function Guesses
    cp = pf0(1,icol);      %Consumption current period
    pigap = pf0(2,icol);   %Inflation gap current period  
    n = pf0(3,icol);       %Labor current period 
    q = pf0(4,icol);       %Tobin's q current period
    %ups = pf0(5,icol);     %Utilization current period
    %----------------------------------------------------------------------
    % Current period
    %----------------------------------------------------------------------
    % HH FOC utilization (1)
    %rk = S.rk*exp(P.sigups*(ups-1));
    % Production function (2)
    yf = (k/g)^P.alpha*n^(1-P.alpha); 
    % Utilization definition (3)
    %u = S.rk*(exp(P.sigups*(ups-1))-1)/P.sigups;
    % Firm FOC capital (4)
    mc = rk*k/(P.alpha*g*yf);
    % Firm FOC labor (5)
    wp = (1-P.alpha)*mc*yf/n;
    % Real wage growth gap (6)
    %wg = pigap*g*wp/(P.g*w);
    % Output definition (7)
    yp = (1-P.varphi*(pigap-1)^2/2)*yf;
    % Lagged ARC (13)
    y = c + x;
    % Output growth gap (8)
    yg = g*yp/(P.g*y);
    % Notional Interest Rate (9)
    inp = in^P.rhoi*(S.i*pigap^P.phipi*yg^P.phiy)^(1-P.rhoi)*exp(mp); 
    % Nominal Interest Rate (10)
%     i = max(1,inp);
    i = inp;
    % Inverse MUC (11)
    lam = cp-P.h*c/g;
    % Flexible real wage definition (12)
    wf = S.chi*n^P.eta*lam;
    % ARC (13)
    xp = yp-cp;  
    % Investment growth gap (14)
    xg = g*xp/(P.g*x);
    % Law of motion for capital (15)
    kp = (1-P.delta)*(k/g)+xp*(1-(xg-1)^2/2);         
    %----------------------------------------------------------------------
    % Linear interpolation of the policy functions 
    %---------------------------------------------------------------------- 
    [cppArr3,pigappArr3,npArr3,qpArr3] = Fallterp743_R(...
        O.g_pts,O.s_pts,O.mp_pts,...
        O.in_pts,O.c_pts,O.k_pts,O.x_pts,O.w_pts,...
        G.in_grid,G.c_grid,G.k_grid,G.x_grid,...
        inp,cp,kp,xp,...
        pf.c,pf.pigap,pf.n,pf.q);
    %----------------------------------------------------------------------        
    % Next period
    %----------------------------------------------------------------------  
    % HH FOC Utilization (1)
    %rkpArr3 = S.rk*exp(P.sigups*(upspArr3-1));  
    % Production function (2)
    yfpArr3 = (kp./gpArr3).^P.alpha.*npArr3.^(1-P.alpha);
    % Utilization definition (3)
    %upArr3 = S.rk*(exp(P.sigups*(upspArr3-1))-1)/P.sigups;
    % Firm FOC capital (4)
    mcpArr3 = rkpArr3.*kp./(P.alpha*gpArr3.*yfpArr3);
    % Firm FOC labor (5)
    wppArr3 = (1-P.alpha)*mcpArr3.*yfpArr3./npArr3;
    % Real wage growth gap (6)
    wgpArr3 = pigappArr3.*gpArr3.*wppArr3/(P.g*wp);
    % Output definition (7)
    yppArr3 = (1-P.varphip*(pigappArr3-1).^2/2).*yfpArr3;
    % Inverse MUC (11)
    lampArr3 = cppArr3-P.h*cp./gpArr3;
    % ARC (13)
    xppArr3 = yppArr3-cppArr3;
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
    %EwpcArr3 = weightArr3.*sdfArr3.*(wgpArr3-1).*wgpArr3.*(yfpArr3/yf);
    Ebond = sum(EbondArr3(:));
    Ecap = sum(EcapArr3(:));
    Einv = sum(EinvArr3(:));
    Eppc = sum(EppcArr3(:));
    %Ewpc = sum(EwpcArr3(:));
    %----------------------------------------------------------------------
    % Euler Equations
    %----------------------------------------------------------------------
    % HH FOC bond (16)
    Res(1,icol) = 1-s*i*Ebond;
    % HH FOC capital (17)
    Res(2,icol) = q-Ecap;
    % HH FOC investment (18)
    Res(3,icol) = 1-q*(1-(xg-1)^2/2-(xg-1)*xg)-P.g*Einv;    
    % Price Phillips Curve (19)
    Res(4,icol) = P.varphi*(pigap-1)*pigap-(1-P.theta)-P.theta*mc-P.varphi*Eppc;
    % Wage Phillips Curve (20)
    %Res(5,icol) = P.varphiw*(wg-1)*wg-((1-P.thetaw)*wp+P.thetaw*wf)*n/yf-P.varphiw*Ewpc;
end
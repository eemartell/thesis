function Res = eqm(pf0,state,O,P,S,G,pf,gpArr3,weightArr3,varargin)

% Get original grids and GH nodes
if ~isempty(varargin)
    EEflag = 1;
    GH = varargin{1};
else
    EEflag = 0;
end

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
    %cp = pf0(1,icol);      %Consumption current period
    pigap = pf0(1,icol);   %Inflation gap current period  
    n = pf0(2,icol);       %Labor current period 
    q = pf0(3,icol);       %Tobin's q current period
    %ups = pf0(5,icol);     %Utilization current period
    mc = pf0(4,icol);      %Marginal cost current period    
    %----------------------------------------------------------------------
    % Current period
    %----------------------------------------------------------------------
    % Production function (2)
    y = (k/g)^P.alpha*n^(1-P.alpha); 
    % Real gdp
    rgdp = c + x;
    rgdpp = (1-P.varphi*(pigap-1)^2/2)*y;
    % Output growth
    rgdpg = g*rgdpp/(P.g*rgdp);    
    % Notional Interest Rate (9)
    inp = in^P.rhoi*(S.i*pigap^P.phipi*rgdpg^P.phiy)^(1-P.rhoi)*exp(mp); 
    % Nominal Interest Rate (10)
    i = inp;    
    % Firm FOC labor (5)
    w = (1-P.alpha)*mc*y/n;
    % FOC labor
    cp = w/(S.chi*n^P.eta)+P.h*c/g;
    % Aggregate resource constraint
    xp = rgdpp - cp;
    % Investment growth gap (14)
    xg = g*xp/(P.g*x);    
    % Law of motion for capital (15)
    kp = (1-P.delta)*(k/g)+xp*(1-P.nu*(xg-1)^2/2);       
    % Inverse MUC (11)
    lam = cp-P.h*c/g;
    %----------------------------------------------------------------------
    % Linear interpolation of the policy functions 
    %---------------------------------------------------------------------- 
    if ~EEflag
    [pigappArr3,npArr3,qpArr3,mcpArr3] = Fallterp743_R(...
        O.g_pts,O.s_pts,O.mp_pts,...
        O.in_pts,O.c_pts,O.k_pts,O.x_pts,...
        G.in_grid,G.c_grid,G.k_grid,G.x_grid,...
        inp,cp,kp,xp,...
        pf.pigap,pf.n,pf.q,pf.mc);
    else
        % Growth rate (9)
        gpVec = P.g + GH.e_nodes;
        % Risk premium (10)
        spVec = (1-P.rhos)*P.s + P.rhos*s + GH.u_nodes;
        [pigappArr3,npArr3,qpArr3,mcpArr3] = Fallterp743(...
            G.g_grid,G.s_grid,G.mp_grid,G.in_grid,...
            G.c_grid,G.k_grid,G.x_grid,...
            inp,cp,kp,xp,...
            gpVec,spVec,GH.v_nodes,...
            pf.pigap,pf.n,pf.q,pf.mc);
        gpArr3 = gpVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1));    
    end
    %----------------------------------------------------------------------        
    % Next period
    %----------------------------------------------------------------------  
    % Production function (2)
    ypArr3 = (kp./gpArr3).^P.alpha.*npArr3.^(1-P.alpha);
    % Firm FOC capital (4)
    rkpArr3 = P.alpha.*mcpArr3.*gpArr3.*ypArr3/kp;
    % Firm FOC labor (5)
    wpArr3 = (1-P.alpha)*mcpArr3.*ypArr3./npArr3;
    % FOC labor
    cppArr3 = wpArr3./(S.chi*npArr3.^P.eta)+P.h*cp./gpArr3;
    % Output definition (7)
    yppArr3 = (1-P.varphi*(pigappArr3-1).^2/2).*ypArr3;
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
    EppcArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*(ypArr3/y);
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
    Res(3,icol) = 1-q*(1-P.nu*(xg-1)^2/2-P.nu*(xg-1)*xg)-P.nu*P.g*Einv;    
    % Price Phillips Curve (19)
    Res(4,icol) = P.varphi*(pigap-1)*pigap-(1-P.theta)-P.theta*mc-P.varphi*Eppc;
    % Wage Phillips Curve (20)
    %Res(5,icol) = P.varphiw*(wg-1)*wg-((1-P.thetaw)*wp+P.thetaw*wf)*n/yf-P.varphiw*Ewpc;
end
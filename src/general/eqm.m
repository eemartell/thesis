function Res = eqm(x,state,O,P,S,G,pf,gpArr3,weightArr3)

% Preallocate function output
Res = zeros(size(x));
ncol = size(Res,2);

% State Values
g = state(1);       %Growth state current period
s = state(2);       %Risk premium state current period
mp = state(3);      %Monetary policy shock current period
in = state(4);  	%Notional interest rate last period
c = state(5);   	%Consumption last period

for icol = 1:ncol   
    % Policy Function Guesses
    cp = x(1,icol);     %Consumption current period 
    pigap = x(2,icol);	%Inflation gap current period     
    %----------------------------------------------------------------------
    % Current period
    %----------------------------------------------------------------------
    % ARC(1) and Output definition (2)
    yf = cp/(1-P.varphip*(pigap-1)^2/2);
    % ARC (1) and Output growth gap (3)
    yg = g*cp/(P.g*c);
    % Nominal interest rate (4)
    inp = in^P.rhoi*(S.i*pigap^P.phipi*yg^P.phiy)^(1-P.rhoi)*exp(mp);
    % Notional interest rate (5)
    i = max(1,inp);
%     i = inp;
    % Inverse MUC (6)
    lam = cp-P.h*c/g;
    % Production function (7) and HH FOC Labor (8)
    w = S.chi*yf^P.eta*lam;
    %----------------------------------------------------------------------
    % Linear interpolation of the policy functions 
    %---------------------------------------------------------------------- 
    [cppArr3,pigappArr3] = Fallterp523_R(...
        O.g_pts,O.s_pts,O.mp_pts,O.in_pts,O.c_pts,...
        G.in_grid,G.c_grid,...
        inp,cp,...
        pf.c,pf.pigap);
    %----------------------------------------------------------------------        
    % Next period
    %----------------------------------------------------------------------    
    % ARC(1) and Output definition (6)
    yfpArr3 = cppArr3./(1-P.varphip*(pigappArr3-1).^2/2);
    % Inverse MUC (5)
    lampArr3 = cppArr3-P.h*cp./gpArr3;
    % Stochastic discount factor
    sdfArr3 = P.beta*lam./lampArr3;
    %----------------------------------------------------------------------
    % Expectations
    %----------------------------------------------------------------------
    EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*(P.pi*pigappArr3));
    EpcArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*(yfpArr3/yf);
    Ebond = sum(EbondArr3(:));
    Epc = sum(EpcArr3(:));
    %----------------------------------------------------------------------
    % Euler Equations
    %----------------------------------------------------------------------
    % HH FOC bond (9)
    Res(1,icol) = 1-s*i*Ebond;
    % Price Phillips Curve (10)
    Res(2,icol) = P.varphip*(pigap-1)*pigap-(1-P.thetap)-P.thetap*w-P.varphip*Epc;
end
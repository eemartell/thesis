function [T,M,eu] = linmodel(P,S,V)

% [T,M,eu] = linmodel(P,S,V)
%   Solves the log-linear model with GENSYS
% Inputs:
%     P     :   Structure of parameters
%     S     :   Structure of steady state values
%     V     :   Structure of variable locations and names
% Output:
%     T     :   Transition matrix
%     M     :   Impact matrix
%     eu    :   [existence, uniqueness]

%---------------------------------------------
%   Initialize GENSYS components
%---------------------------------------------
G0 = zeros(V.nvar);
G1 = zeros(V.nvar);
Psi = zeros(V.nvar,V.nshock);
Pi = zeros(V.nvar,V.nfore);
CC = zeros(V.nvar,1);
j = 0;
%---------------------------------------------
j=j+1;%	[FOC Labor]
%---------------------------------------------
G0(j,V.w) = 1;
G0(j,V.n) = -P.eta;
G0(j,V.c) = -1/(1-P.htilde);
G0(j,V.g) = -P.htilde/(1-P.htilde);
G1(j,V.c) = -P.htilde/(1-P.htilde);
%---------------------------------------------
j = j+1;% [FOC Bond]
%---------------------------------------------
G0(j,V.c) = 1/(1-P.htilde);
G0(j,V.pi) = 1;
G0(j,V.beta) = -1;
G0(j,V.g) = 1/(1-P.htilde);
G1(j,V.r) = 1;
G1(j,V.c) = (1+P.htilde)/(1-P.htilde);
G1(j,V.g) = P.htilde/(1-P.htilde);
G1(j,V.clag) = -P.htilde/(1-P.htilde);
Pi(j,V.fec) = 1/(1-P.htilde);
Pi(j,V.fepi) = 1;
Pi(j,V.febeta) = -1;
Pi(j,V.feg) = 1/(1-P.htilde);
%---------------------------------------------
j = j+1;% [FOC Capital]
%---------------------------------------------
G0(j,V.q) = P.beta*(1-P.delta)/P.g;
G0(j,V.beta) = 1;
G0(j,V.rk) = 1-P.beta*(1-P.delta)/P.g;
G0(j,V.g) = -1/(1-P.htilde);
G0(j,V.c) = -1/(1-P.htilde);
G1(j,V.q) = 1;
G1(j,V.c) = -(1+P.htilde)/(1-P.htilde);
G1(j,V.g) = -P.htilde/(1-P.htilde);
G1(j,V.clag) = P.htilde/(1-P.htilde);
Pi(j,V.feq) = P.beta*(1-P.delta)/P.g;
Pi(j,V.febeta) = 1;
Pi(j,V.ferk) = 1-P.beta*(1-P.delta)/P.g;
Pi(j,V.feg) = -1/(1-P.htilde);
Pi(j,V.fec) = -1/(1-P.htilde);
%---------------------------------------------
j = j+1;% [FOC Investment]
%---------------------------------------------
G0(j,V.i) = P.nu*P.beta;
G0(j,V.g) = P.nu*P.beta;
G1(j,V.q) = -1;        
G1(j,V.g) = P.nu;
G1(j,V.i) = P.nu*(1+P.beta);
G1(j,V.ilag) = -P.nu;
Pi(j,V.fei) = P.nu*P.beta;
Pi(j,V.feg) = P.nu*P.beta;
%---------------------------------------------
j = j+1;% [Production Function] 
%---------------------------------------------
G0(j,V.y) = 1;
G0(j,V.n) = -(1-P.alpha);
G0(j,V.g) = P.alpha;
G1(j,V.k) = P.alpha;
%---------------------------------------------
j = j+1;% [Law of Motion for Capital]
%---------------------------------------------
G0(j,V.k) = S.k;
G0(j,V.i) = -S.i;
G0(j,V.g) = (1-P.delta)*S.k/P.g;
G1(j,V.k) = (1-P.delta)*S.k/P.g;
%---------------------------------------------
j = j+1;% [Firm FOC n] 
%---------------------------------------------
G0(j,V.w) = 1;
G0(j,V.y) = -1;
G0(j,V.n) = 1;
G0(j,V.psi) = -1;
%---------------------------------------------
j = j+1;% [Firm FOC k] 
%---------------------------------------------
G0(j,V.rk) = 1;
G0(j,V.y) = -1;
G0(j,V.g) = -1;
G0(j,V.psi) = -1;
G1(j,V.k) = -1;
%---------------------------------------------
j = j+1;% [Firm Pricing] 
%---------------------------------------------
G0(j,V.pi) = P.beta;
G1(j,V.pi) = 1;
G1(j,V.psi) = (1-P.theta)/P.varphi;
Pi(j,V.fepi) = P.beta; 
%---------------------------------------------
j = j+1;% [ARC]
%---------------------------------------------
G0(j,V.c) = S.c;
G0(j,V.i) = S.i;
G0(j,V.y) = -S.y;
%---------------------------------------------
j=j+1;  %   [Notional Interest Rate] 
%---------------------------------------------
G0(j,V.rn) = 1;
G0(j,V.r) = -1;
%---------------------------------------------
j = j+1;% [Interest Rate Rule]
%---------------------------------------------
G0(j,V.rn) = 1;
G0(j,V.pi) = -P.phipi*(1-P.rhor);
G0(j,V.rgdp) = -P.phiy*(1-P.rhor);
G0(j,V.g) = -P.phiy*(1-P.rhor);
G1(j,V.rn) = P.rhor;
G1(j,V.rgdp) = -P.phiy*(1-P.rhor);
Psi(j,V.epsmp) = 1;
%---------------------------------------------
j=j+1;  %   [Growth Process] 
%---------------------------------------------
G0(j,V.g) = 1;
G1(j,V.g) = P.rhog;
Psi(j,V.epsg) = 1;
%---------------------------------------------
j=j+1;  %   [Discount Factor Shock Process] 
%---------------------------------------------
G0(j,V.beta) = 1;
G1(j,V.beta) = P.rhobeta;
Psi(j,V.epsbeta) = 1;
%---------------------------------------------
j=j+1;  %   [Real GDP] 
%---------------------------------------------
G0(j,V.c) = S.c;
G0(j,V.i) = S.i;
G0(j,V.rgdp) = -S.rgdp;
%---------------------------------------------
j=j+1;  %   [Lagged Consumption] 
%---------------------------------------------
G0(j,V.clag) = 1;
G1(j,V.c) = 1;
%---------------------------------------------
j=j+1;  %   [Lagged Investment] 
%---------------------------------------------
G0(j,V.ilag) = 1;
G1(j,V.i) = 1;

%---------------------------------------------
%   Solve Linear Model
%---------------------------------------------
[T,~,M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);        
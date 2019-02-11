function [T, M, eu] = linmodel(P,S,V)

% [T, M, eu] = linmodel(P,S,V)
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
G0(j,V.c) = -P.sigma;
G0(j,V.n) = -P.eta;
%---------------------------------------------
j = j+1;% [Consumption Euler Equation]
%---------------------------------------------
G0(j,V.c) = P.sigma;
G0(j,V.rk) = -(1-P.beta*(1-P.delta));
G0(j,V.beta) = -1;
G0(j,V.q) = -P.beta*(1-P.delta);
G0(j,V.i) = -P.beta*P.nu*P.delta^2;
G1(j,V.q) = -1;
G1(j,V.c) = P.sigma;
G1(j,V.k) = -P.beta*P.nu*P.delta^2;
Pi(j,V.fec) = P.sigma;
Pi(j,V.ferk) = -(1-P.beta*(1-P.delta));
Pi(j,V.febeta) = -1;
Pi(j,V.feq) = -P.beta*(1-P.delta);
Pi(j,V.fei) = -P.beta*P.nu*P.delta^2;
%---------------------------------------------
j = j+1;% [FOC Bond]
%---------------------------------------------
G0(j,V.c) = P.sigma;
G0(j,V.pi) = 1;
G0(j,V.beta) = -1;
G1(j,V.r) = 1;
G1(j,V.c) = P.sigma;
Pi(j,V.fec) = P.sigma;
Pi(j,V.fepi) = 1;
Pi(j,V.febeta) = -1;
%---------------------------------------------
j = j+1;% [ARC]
%---------------------------------------------
G0(j,V.c) = S.c;
G0(j,V.y) = -S.y;
G0(j,V.i) = S.i;
%---------------------------------------------
j = j+1;% [Production Function] 
%---------------------------------------------
G0(j,V.y) = 1;
G0(j,V.n) = -(1-P.alpha);
G1(j,V.k) = P.alpha;
G0(j,V.z) = -1;
%---------------------------------------------
j = j+1;% [Firm Pricing] 
%---------------------------------------------
G0(j,V.pi) = P.beta;
G1(j,V.pi) = 1;
G1(j,V.psi) = (1-P.theta)/P.varphi;
Pi(j,V.fepi) = P.beta; 
%---------------------------------------------
j = j+1;% [Marginal Cost Definition]
%---------------------------------------------
G0(j,V.psi) = 1;
G0(j,V.w) = -(1-P.alpha);
G0(j,V.rk) = -P.alpha;
G0(j,V.z) = 1;
%---------------------------------------------
j = j+1;% [Interest Rate Rule]
%---------------------------------------------
G0(j,V.r) = 1;
G0(j,V.pi) = -P.phipi;
G0(j,V.y) = -P.phiy;
%---------------------------------------------
j = j+1;% [Law of Motion for Capital]
%---------------------------------------------
G0(j,V.i) = P.delta;
G0(j,V.k) = -1;
G1(j,V.k) = -(1-P.delta);
%---------------------------------------------
j = j+1;% [Firm Pricing] 
%---------------------------------------------
G0(j,V.w) = 1; 
G0(j,V.n) = 1;
G0(j,V.rk) = -1;
G1(j,V.k) = 1;
%---------------------------------------------
j=j+1;  %   [Productivity Process] 
%---------------------------------------------
G0(j,V.z) = 1;
G1(j,V.z) = P.rhoz;
Psi(j,V.epsz) = 1;
%---------------------------------------------
j=j+1;  %   [Tobin's q] 
%---------------------------------------------
G0(j,V.q) = 1;
G0(j,V.i) = -P.nu*P.delta;
G1(j,V.k) = -P.nu*P.delta;
%---------------------------------------------
j=j+1;  %   [Discount factor Process] 
%---------------------------------------------
G0(j,V.beta) = 1;
G1(j,V.beta) = P.rhobeta;
Psi(j,V.epsbeta) = 1;

%---------------------------------------------
%   Solve Linear Model
%---------------------------------------------
[T,~,M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);        
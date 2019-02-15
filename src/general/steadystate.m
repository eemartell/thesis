function S = steadystate(P)

% S = steadystate(P) 
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% FOC bond 
S.i = P.g*P.pi/(P.beta*P.s);
S.in = S.i;
% Firm pricing
S.mc = (P.theta-1)/P.theta;
% FOC capital
S.rk = P.g/P.beta+P.delta-1;
% Marginal cost definition
S.w = (S.mc*(1-P.alpha)^(1-P.alpha)*P.alpha^P.alpha/S.rk^P.alpha)^(1/(1-P.alpha));
%S.wf = (P.thetaw-1)*S.w/P.thetaw;
% Consolidated FOC firm
S.k = S.w*P.n*P.g*P.alpha/(S.rk*(1-P.alpha));
% Law of motion for capital
S.x = (1-(1-P.delta)/P.g)*S.k;
% Production function
S.yf = (S.k/P.g)^P.alpha*P.n^(1-P.alpha);
% Real GDP
S.y = S.yf;
% Aggregate resouce constraint
S.c = S.y-S.x;
% FOC labor
S.lam = (1-P.h/P.g)*S.c;
S.chi = S.w/(P.n^P.eta*S.lam);
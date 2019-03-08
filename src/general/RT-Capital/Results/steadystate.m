function S = steadystate(P)

% S = steadystate(P) 
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% FOC bond 
S.r = P.pi*P.g/P.beta;
% Firm pricing
S.psi = 1/P.mu;
% FOC investment
S.q = 1;
% FOC capital
S.rk = S.q*(P.g/P.beta+P.delta-1);
% Marginal cost definition
S.w = (S.psi*(1-P.alpha)^(1-P.alpha)*P.alpha^P.alpha/S.rk^P.alpha)^(1/(1-P.alpha));
% Consolidated FOC firm
S.k = S.w*P.n*P.g*P.alpha/(S.rk*(1-P.alpha));
% Law of motion for invesatment
S.i = (1-(1-P.delta)/P.g)*S.k;
% Production function
S.y = (S.k/P.g)^P.alpha*P.n^(1-P.alpha);
% Aggregate resouce constraint
S.c = S.y-S.i;
% Real GDP
S.rgdp = S.c+S.i;
% FOC labor
S.chi = S.w/(P.n^P.eta*S.c);
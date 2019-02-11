function S = steadystate(P)

% S = steadystate(P)
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% FOC bond 
S.r = P.pi/P.beta;
% FOC capital
S.rk = 1/P.beta + P.delta - 1;
% Firm pricing
S.psi = (P.theta-1)/P.theta;
% Marginal cost definition
S.w = (S.psi*(1-P.alpha)^(1-P.alpha)*P.alpha^P.alpha/S.rk^P.alpha)^(1/(1-P.alpha));
% Consolidated FOC firm
S.k = S.w*P.n*P.alpha/(S.rk*(1-P.alpha));
% Investment definition
S.i = P.delta*S.k;
% Prod. Tech
S.y = P.zbar*S.k^P.alpha*P.n^(1-P.alpha);
% Aggregate resouce constraint
S.c = S.y - S.i;
% FOC labor
S.chi = S.w/(P.n^P.eta*S.c^P.sigma);
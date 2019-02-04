function P = parameters

% P = parameters
%   Defines model parameters
% Output:
%   P : Structure of model parameters

% Algorithm
P.tol      = 1e-6;      % Convergence criterion

% Calibrated Parameters
P.beta     = 0.9955;    % Discount factor
P.theta    = 6;         % Elasticity of substitution 
P.n        = 1/3;       % Steady state labor
P.eta      = 1/3;       % Inverse Frisch elasticity of labor supply
P.rhog     = 0;         % Productivity growth persistence
P.rhomp    = 0;  	    % Monetary policy shock persistence 
P.s        = 0.006;    % Average risk premium
P.g        = 1.0034;    % Average growth rate
P.pi       = 1.0048;    % Inflation target

% Estimated Parameters
P.varphi   = 100;       % Price adjustment cost
P.phipi    = 2;         % Inflation response
P.rhoa     = 0.80;      % Risk premium shock persistence
P.rhoi     = 0.75;      % Notional rate persistence
P.sige     = 0.005;  	% Productivity growth shock std
P.sigu     = 0.005;  	% Preference shock std
P.sigv     = 0.005;  	% Interest rate shock std

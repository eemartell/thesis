function P = parameters

% Calibrated Parameters
P.beta     = 0.9949;    % Discount factor
P.theta    = 6;
P.n        = 1/3;       % Steady state labor
P.eta      = 1/3;       % Inverse frish elasticity of labor supply
P.delta    = 0.025;     % Depreciation
P.alpha    = 0.35;      % Capital share
P.g        = 1.0034;    % Mean growth rate
P.pi       = 1.0053;   	% Inflation target
P.s        = 1.0058;    % Average risk premium

% Parameters for DGP and Estimated parameters
P.varphi  = 100;       % Rotemberg price adjustment cost
P.h        = 0.80;  	% Habit persistence
P.rhos     = 0.80;      % Persistence
P.rhoi     = 0.80; 	    % Persistence
P.sigg     = 0.005;     % Standard deviation
P.sigs     = 0.005;     % Standard deviation
P.sigmp    = 0.002;     % Standard deviation
P.phipi    = 2.0;       % Inflation responsiveness
P.phiy     = 0.5;       % Output responsiveness
P.nu       = 4;         % Investment adjustment cost

% Algorithm
P.tol      = 1e-6;      % Convergence criterion
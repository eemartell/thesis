function P = parameters

% Calibrated Parameters
P.beta     = 0.9983;    % Discount factor
P.theta    = 6;         % Price elasticity of demand
P.n        = 1/3;       % Steady state labor
P.zlb      = 1.00022;	% Lower bound
P.eta      = 1/3;       % Inverse frish elasticity of labor supply
P.delta    = 0.025;     % Depreciation

% Estimated Parameters
P.varphi   = 95;        % Rotemberg parameter
% P.h        = 0.5;  	    % Habit persistence
P.rhobeta  = 0.90;      % Persistence
P.rhog     = 0.1;       % Persistence
P.rhor     = 0.8; 	    % Persistence
P.sige     = 0.007;     % Standard deviation of shock
P.sigu     = 0.0025;    % Standard deviation
P.sigv     = 0.001;     % Standard deviation
P.phipi    = 3.5;       % Inflation responsiveness
P.phiy     = 0.5;       % Output responsiveness
P.g        = 1.0023;    % Mean
P.pi       = 1.0057;   	% Inflation target
P.alpha    = 0.33;      % Capital share
P.nu       = 4;         % Investment adjustment cost

% Implied parameters
P.mu       = P.theta/(P.theta-1);
% P.htilde   = P.h/P.g;

% Algorithm
P.tol      = 1e-6;      % Convergence criterion
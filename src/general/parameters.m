function P = parameters

% P = parameters
%   Defines model parameters
% Output:
%   P : Structure of model parameters

% Household (Quarterly Calibration)
P.beta     = 0.995;    % Discount factor
P.sigma    = 1;        % Constant of relative risk aversion
P.eta      = 1/3;      % Inverse frish elasticity of labor supply
P.theta    = 6;        % Price elasticity of demand
P.n        = 1/3;      % Steady state labor
P.omega    = 0.75;     % Calvo parameter
P.nu       = 5.6;      % Capital adjustment cost coefficient

% Firm
P.alpha    = 0.33;     % Capital share
P.delta    = 0.025;    % Depreciation

% Technology Process
P.zbar     = 1;        % Mean
P.rhoz     = 0.90;     % Persistence
P.sige     = 0.0025;   % Standard deviation

% Discount Factor Process
P.rhobeta  = 0.80;     % Persistence
P.sigu     = 0.0025;   % Standard deviation

% Monetary Policy 
P.pi       = 1.006;    % Inflation target
P.phipi    = 1.5;      % Inflation response
P.phiy     = 0.05;     % Output response

% Implied parameters
P.varphi   = P.omega*(P.theta-1)/((1-P.omega)*(1-P.beta*P.omega));
P.mu       = P.theta/(P.theta-1);

% Algorithm
P.tol      = 1e-10;    % Convergence criterion
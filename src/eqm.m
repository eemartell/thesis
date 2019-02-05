function R = eqm(x,state,O,P,S,G,pf,gpArr3,apArr3,weightArr3,varargin)

% Get original grids and GH nodes
if ~isempty(varargin)
    EEflag = 1;
    GH = varargin{1};
else
    EEflag = 0;
end

% Preallocate function output
R = zeros(size(x));
Rdim = size(R,2);

% State Values
g = state(1);       %Growth state current period
a = state(2);       %Preference state current period
mp = state(3);      %Monetary policy state current period
in = state(4);  	%Notional interest rate last period

for icol = 1:Rdim   
    % Policy Function Guesses
    cp = x(1,icol);     %Consumption policy current period 
    pigap = x(2,icol);	%Inflation policy current period     
    %----------------------------------------------------------------------
    % Solve for variables
    %----------------------------------------------------------------------
    % Aggregate resource constraint (6)    
    y = cp/(1-P.varphi*(pigap-1)^2/2);
    % Interest rate rule (8)
    inp = in^P.rhoi*(S.i*pigap^P.phipi)^(1-P.rhoi)*exp(mp);
    
    if P.zlbflag
        i = max(1,inp);
    else
        i = inp;
    end
    % FOC Labor (1,2)
    lam = cp/a;
    w = S.chi*a*y^P.eta*lam;
    %----------------------------------------------------------------------
    % Linear interpolation of the policy variables
    %----------------------------------------------------------------------
    if ~EEflag 
	    [cppArr3,pigappArr3] = Fallterp423_R(...
		O.g_pts,O.a_pts,O.mp_pts,O.in_pts,...
		G.in_grid,...
		inp,...
		pf.c,pf.pigap);
    else
        % Growth rate
        gpVec = ((1-P.rhog)*P.g + P.rhog*g + GH.e_nodes);
        % Preferences 
        apVec = (1-P.rhoa + P.rhoa*a + GH.u_nodes);
        [cppArr3,pigappArr3] = allterp423(...
                                G.g_grid,G.a_grid,G.mp_grid,G.in_grid,...
                                gpVec,apVec,GH.v_nodes,inp,...
                                pf.c,pf.pigap);                            
        gpArr3 = gpVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1));
        apArr3 = permute(apVec(:,ones(GH.shockpts,1),ones(GH.shockpts,1)),[2,1,3]);
    end
    %----------------------------------------------------------------------        
    % Solve for variables inside expectations
    %----------------------------------------------------------------------    
    % Aggregate resource constraint  
    ypArr3 = cppArr3./(1-P.varphi*(pigappArr3-1).^2/2);
    % Stochastic discount factor
    lampArr3 = cppArr3./apArr3;
    sdfArr3 = P.beta*lam./lampArr3;
    %----------------------------------------------------------------------
    % Numerical integration
    %----------------------------------------------------------------------
    % Compute all combinations of shocks
    EbondArr3 = weightArr3.*sdfArr3./(gpArr3.*pigappArr3);
    EfpArr3 = weightArr3.*sdfArr3.*(pigappArr3-1).*pigappArr3.*ypArr3;
    % Integrate
    Ebond = sum(EbondArr3(:));
    Efp = sum(EfpArr3(:));
    %----------------------------------------------------------------------
    % First-order conditions
    %----------------------------------------------------------------------
    % Consumption Euler equation (3)
    R(1,icol) = 1 - (1+P.s)*i*Ebond/P.pi;
    % Firm Pricing (Philips Curve, 5)
    R(2,icol) = P.varphi*(pigap-1)*pigap-(1-P.theta)-P.theta*w-P.varphi*Efp/y;
end

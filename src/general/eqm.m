function R = eqm(x,state,...
                Ozlbflag,...
                Palpha,Ppi,Pphipi,Pphiy,Peta,Pvarphi,...  
                Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,Ptheta,...
                Schi,Sr,Sy,...       
                Ge_weight,Ge_nodes,Gu_weight,Gu_nodes,...
                Gk_grid,Gz_grid,Gbeta_grid,...
                pf_n,pf_pi,pf_i)

% R = eqm(x,state,P,S,G,pf,zlbflag)
%   Outputs residuals of the equilibrium system of equations for time
%   iteration/linear interpolation method.
% Inputs:
%   x       :   Policy function values on node i
%   state   :   State variable value on node i
%   O       :   Structure of options
%   P       :   Structure of parameters
%   S       :   Structure of steady state values
%   G       :   Structure of grids
%   pf      :   Structure of policy functions
% zlbflag   :   Flag to impose zero-lower bound on the interest rate
% Output:
%   R       :   Vector of residuals on node i

% Preallocate function output
R = zeros(size(x));
Rdim = size(R,2);

% State Values
k = state(1);       %Capital state last period  
z = state(2);       %Technology state current period
beta = state(3);    %Discount factor state current period

for j = 1:Rdim   
    % Policy Function Guesses
    n = x(1,j);     %Labor policy current period 
    pie = x(2,j);   %Inflation policy current period 
    i = x(3,j); 	%Investment policy current period 
    %----------------------------------------------------------------------
    % Solve for variables
    %----------------------------------------------------------------------
    % Production function
    y = z*k^Palpha*n^(1-Palpha);    
    % Interest rate rule
    r = Sr*(pie/Ppi)^Pphipi*(y/Sy)^Pphiy;   
    if Ozlbflag == 1
        r = max(1,r);
    end       
    % Tobin's q
    q = 1+Pnu*(i/k-Pdelta);    
    % Aggregate resource constraint
    ytil = y*(1-(Pvarphi*(pie/Ppi-1)^2)/2);
    kac = Pnu*(i/k-Pdelta)^2/2;
    c = ytil - i - kac*k;    
    % FOC Labor
    w = Schi*n^Peta*c;
    % Firm FOC labor
    psi = w*n/((1-Palpha)*y);
    % Investment
    kp = i + (1-Pdelta)*k;       
    % Technology
    zpVec = Pzbar*(z/Pzbar)^Prhoz*exp(Ge_nodes); 
    % Discount factor
    betapVec = Pbeta*(beta/Pbeta)^Prhobeta*exp(Gu_nodes);      
    %----------------------------------------------------------------------
    % Linear interpolation of the policy variables 
    %---------------------------------------------------------------------- 
    [npMat,pipMat,ipMat] = Fallterp332c(Gk_grid,Gz_grid,Gbeta_grid,...
                                        kp,zpVec,betapVec,...
                                        pf_n,pf_pi,pf_i);    
    %----------------------------------------------------------------------        
    % Solve for variables inside expectations
    %----------------------------------------------------------------------           
    zpMat = zpVec(:,ones(numel(Gu_weight),1));
    betapMat = betapVec(:,ones(length(Ge_nodes),1))';    
    ypMat = zpMat.*(kp^Palpha).*(npMat.^(1-Palpha));      
    qpMat = 1+Pnu*(ipMat/kp-Pdelta);  
    ytilpMat = ypMat.*(1-(Pvarphi*(pipMat/Ppi-1).^2)/2);
    kacpMat = Pnu*(ipMat/kp-Pdelta).^2/2;
    cpMat = ytilpMat - ipMat - kacpMat*kp;      
    wpMat = Schi*npMat.^Peta.*cpMat;   
    psipMat = wpMat.*npMat./((1-Palpha)*ypMat);
    rkpMat = Palpha*psipMat.*ypMat/kp;
    sdfMat = betapMat.*(c./cpMat);   
    %----------------------------------------------------------------------
    % Numerical integration
    %----------------------------------------------------------------------    
    eMat = Ge_weight(:,ones(numel(Gu_weight),1));
    % Apply GH across e
    Efp =  sum(eMat.*(sdfMat.*(pipMat./Ppi-1).*(ypMat/y).*(pipMat./Ppi)))/sqrt(pi);
    Ebond = sum(eMat.*(r*sdfMat./pipMat))/sqrt(pi);
    Econs = sum(eMat.*(sdfMat.*(rkpMat-kacpMat+(qpMat-1).*(ipMat/kp)+(1-Pdelta)*qpMat)))/sqrt(pi);
    % Apply GH across u
    Efp = sum(Gu_weight'.*Efp)/sqrt(pi);
    Ebond = sum(Gu_weight'.*Ebond)/sqrt(pi);
    Econs = sum(Gu_weight'.*Econs)/sqrt(pi);     
    %----------------------------------------------------------------------
    % First-order conditions
    %----------------------------------------------------------------------
    % Firm Pricing Equation
    R(1,j) = Pvarphi*(pie/Ppi-1)*pie/Ppi-(1-Ptheta)-Ptheta*psi-Pvarphi*Efp;
    % Bond Euler Equation
    R(2,j) = 1-Ebond;
    % Consumption Euler Equation
    R(3,j) = q-Econs;
end
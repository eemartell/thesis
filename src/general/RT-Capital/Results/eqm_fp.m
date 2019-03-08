function x_up = eqm_fp(x,state,P,S,G,pf,...
                 gpArr,betapArr,weightArr)

% State Values
c = state(1);   	%Consumption last period
k = state(2);  	    %Capital last period
i = state(3);       %Investment last period   
rn = state(4);  	%Notional interest rate last period
g = state(5);       %Growth state current period
mp = state(6);      %Interest rate shock current period 

% Policy Function Guesses
n = x(1);     %Labor policy current period 
q = x(2);     %Tobin's q current period
pie = x(3);   %Inflation current period  
psi = x(4);   %Marginal cost current period
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
% Production function
y = (k/g)^P.alpha*n^(1-P.alpha);
% Real GDP
rgdp = c+i;
rgdpp = (1-(P.varphi*(pie/P.pi-1)^2)/2)*y;
% Interest rate rule
rnp = rn^P.rhor*(S.r*(pie/P.pi)^P.phipi*...
      (g*rgdpp/(P.g*rgdp))^P.phiy)^(1-P.rhor)*exp(mp);   
r = max(P.zlb,rnp);
% Firm FOC labor
w = (1-P.alpha)*psi*y/n;
% FOC labor
cp = w/(S.chi*n^P.eta)+P.h*c/g;
% Aggregate resource constraint
ip = rgdpp - cp;  
% Investment adjustment costs
iac = ip*g/(i*P.g);
% Law of motion for capital
kp = (1-P.delta)*(k/g)+ip*(1-P.nu*(iac-1)^2/2);         
%----------------------------------------------------------------------
% Linear interpolation of the policy variables 
%---------------------------------------------------------------------- 
[npArr,qpArr,pipArr,psipArr] = Fallterp743(...
   G.c_grid,G.k_grid,G.i_grid,G.rn_grid,...
   G.g_grid,G.beta_grid,G.mp_grid,...
   cp,kp,ip,rnp,G.e_nodes,G.u_nodes,G.v_nodes,...
   pf.n,pf.q,pf.pi,pf.psi);
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
% Production function
ypArr = (kp./gpArr).^P.alpha.*npArr.^(1-P.alpha);
% Firm FOC capital
rkpArr = P.alpha*psipArr.*gpArr.*ypArr/kp;
% Firm FOC labor
wpArr = (1-P.alpha)*psipArr.*ypArr./npArr;
% FOC labor
cppArr = wpArr./(S.chi*npArr.^P.eta)+P.h*cp./gpArr;
% Aggregate resource constraint  
ippArr = ypArr.*(1-(P.varphi*(pipArr/P.pi-1).^2)/2) - cppArr;       
% Investment adjustment costs
iacpArr = ippArr.*gpArr/(ip*P.g);
% Stochastic discount factor
sdfArr = betapArr.*(cp-P.h*c/g)./(cppArr-P.h*cp./gpArr);
%----------------------------------------------------------------------
% Numerical integration
%----------------------------------------------------------------------
% Compute all combinations of shocks
EcapArr = weightArr.*(sdfArr./gpArr).*(rkpArr+(1-P.delta)*qpArr);
EinvArr = weightArr.*sdfArr.*qpArr.*(gpArr/P.g).*(ippArr/ip).^2.*(iacpArr-1);
EbondArr = weightArr.*sdfArr./(gpArr.*pipArr);
EfpArr = weightArr.*sdfArr.*(pipArr./P.pi-1).*(ypArr/y).*(pipArr./P.pi);
% Integrate
Ecap = sum(EcapArr(:));
Einv = sum(EinvArr(:));    
Ebond = sum(EbondArr(:));
Efp = sum(EfpArr(:));
%----------------------------------------------------------------------
% First-order conditions
%----------------------------------------------------------------------
x_up(2) = Ecap; %q
RHS = (1-P.theta) + P.theta*mc + P.varphi*Efp;
x_up(3) = 1/2*(1-sqrt(4*RHS+1))*P.pi;


% % FOC capital
% R(1,j) = q - Ecap;
% % FOC investment
% R(2,j) = 1-q*(1-P.nu*(iac-1)^2/2-P.nu*iac*(iac-1))-P.nu*Einv;    
% % FOC bond
% R(3,j) = 1 - r*Ebond;
% % Firm Pricing
% R(4,j) = P.varphi*(pie/P.pi-1)*pie/P.pi-(1-P.theta)-P.theta*psi-P.varphi*Efp;
end
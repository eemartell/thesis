function [Einvpi,Ec,Er,Erk] = expvar(k,z,beta,i,...
                                     Palpha,Ppi,Peta,Pvarphi,Pphipi,Pphiy,...  
                                     Pdelta,Pnu,Pzbar,Prhoz,Pbeta,Prhobeta,...
                                     Schi,Sr,Sy,...       
                                     Ge_weight,Ge_nodes,Gu_weight,Gu_nodes,...
                                     Gk_grid,Gz_grid,Gbeta_grid,...
                                     pf_n,pf_pi,pf_i)

% Replicate state across the shocks                                 
nsims = numel(k);
kArr = k(:,ones(length(Ge_nodes),1),ones(length(Gu_nodes),1));
zArr = z(:,ones(length(Ge_nodes),1),ones(length(Gu_nodes),1));
betaArr = beta(:,ones(length(Ge_nodes),1),ones(length(Gu_nodes),1));
iArr = i(:,ones(length(Ge_nodes),1),ones(length(Gu_nodes),1));
eNodesArr = permute(Ge_nodes(:,ones(nsims,1),ones(length(Gu_nodes),1)),[2 1 3]);
uNodesArr = permute(Gu_nodes(:,ones(nsims,1),ones(length(Ge_nodes),1)),[2 3 1]);
npArr = 0*zArr;
pipArr = 0*zArr;
ipArr = 0*zArr;
%----------------------------------------------------------------------
% Solve for variables
%----------------------------------------------------------------------
% Investment
kpArr = iArr + (1-Pdelta)*kArr;   
% Technology
zpArr = Pzbar*(zArr/Pzbar).^Prhoz.*exp(eNodesArr);
% Discount factor
betapArr = Pbeta*(betaArr/Pbeta).^Prhobeta.*exp(uNodesArr); 
%----------------------------------------------------------------------
% Linear interpolation of the policy variables 
%---------------------------------------------------------------------- 
[npArr(:),pipArr(:),ipArr(:)] = Fallterp33c_F(Gk_grid,Gz_grid,Gbeta_grid,...
                                              kpArr(:),zpArr(:),betapArr(:),...
                                              pf_n,pf_pi,pf_i);  
%----------------------------------------------------------------------        
% Solve for variables inside expectations
%----------------------------------------------------------------------    
% Production function
ypArr = zpArr.*(kpArr.^Palpha).*(npArr.^(1-Palpha));
% Adjusted output
ytilp = ypArr.*(1-(Pvarphi*(pipArr/Ppi-1).^2)/2);
% Capital adjustment costs
kacp = Pnu*(ipArr./kpArr-Pdelta).^2/2;
% Aggregate resource constraint
cpArr = ytilp - ipArr - kacp.*kpArr;    
% FOC Labor
wpArr = Schi*npArr.^Peta.*cpArr;   
% Firm FOC labor
psipArr = wpArr.*npArr./((1-Palpha)*ypArr);
% Firm FOC capital
rkpArr = Palpha*psipArr.*ypArr./kpArr;
% Interest rate rule
rpArr = max(1,Sr*(pipArr/Ppi).^Pphipi.*(ypArr./Sy).^Pphiy);  
%----------------------------------------------------------------------
% Numerical integration
%----------------------------------------------------------------------
eWeightArr = permute(Ge_weight(:,ones(nsims,1),ones(length(Gu_nodes),1)),[2 1 3]);
uWeightArr = permute(Gu_weight(:,ones(nsims,1)),[2 1]);
% Expectations
Einvpi =  eWeightArr.*(1./pipArr);
Ec =  eWeightArr.*cpArr;
Er =  eWeightArr.*rpArr;
Erk =  eWeightArr.*rkpArr;
% Apply GH across e
Einvpi = sum(Einvpi,2)/sqrt(pi);
Ec = sum(Ec,2)/sqrt(pi);
Er = sum(Er,2)/sqrt(pi);
Erk = sum(Erk,2)/sqrt(pi);
% Apply GH across u
if nsims == 1
    Einvpi = uWeightArr'.*squeeze(Einvpi);
    Einvpi = sum(Einvpi)/sqrt(pi); 
    Ec = uWeightArr'.*squeeze(Ec);
    Ec = sum(Ec)/sqrt(pi);
    Er = uWeightArr'.*squeeze(Er);
    Er = sum(Er)/sqrt(pi);    
    Erk = uWeightArr'.*squeeze(Erk);
    Erk = sum(Erk)/sqrt(pi);        
else
    Einvpi = uWeightArr.*squeeze(Einvpi);
    Einvpi = sum(Einvpi,2)/sqrt(pi);
    Ec = uWeightArr.*squeeze(Ec);
    Ec = sum(Ec,2)/sqrt(pi);
    Er = uWeightArr.*squeeze(Er);
    Er = sum(Er,2)/sqrt(pi);    
    Erk = uWeightArr.*squeeze(Erk);
    Erk = sum(Erk,2)/sqrt(pi);    
end
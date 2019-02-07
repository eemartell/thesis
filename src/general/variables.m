function V = variables

% Variables: Names variables and assigns locations
% Inputs:
%     model
% Output:
%     V     :   Structure of variable locations and names

% Variables names
V.names =      {'c'          %1  Consumption 
                'n'          %2  Labor
                'y'          %3  Output
                'yf'         %4  Final goods
                'yg'         %5  Output growth gap
                'w'          %6  Real wage rate
                'pi'         %7  Inflation rate
                'i'          %8  Nominal interest rate
                'in'         %9  Notional interest rate
                'lam'        %10 Inverse MUC
                'g'          %11 Growth shock       
                's'          %12 Risk premium shock
                'mp'         %13 Monetary policy shock    
               };

% Variables titles        
V.desc =      { 'Consumption'
                'Labor'
                'Output'
                'Final Goods'
                'Output Growth Gap'
                'Real Wage'
                'Inflation Rate'
                'Nominal Interest Rate'
                'Notional Interest Rate'
                'Inverse MUC'
                'Growth Shock'
                'Risk Premium Shock'
                'Monetary Policy Shock'                          
            };

% Variables names for nonlinear model
V.plotnames = V.names;
V.nplotvar = length(V.plotnames);
      
% Shocks
V.shocktypes = {'g','s','mp'};
V.nshock = length(V.shocktypes);

% Forecast errors
V.foretypes = {'lam','g','pi'}; 
V.nfore = length(V.foretypes);
             
% Add expectation variables
for j = 1:V.nfore
    V.names = [V.names; ['e' V.foretypes{j}]];
end
V.nvar = length(V.names);

% Establish variable index
%  Linear model
for j = 1:V.nvar
   eval(['V.' V.names{j} ' = j;']);      
end
%  Nonlinear model
for j = V.nvar-V.nfore+1:V.nplotvar
   eval(['V.' V.plotnames{j} ' = j;']);      
end
% Establish shock index
for j = 1:V.nshock
    eval(['V.eps' V.shocktypes{j} ' = j;'])
end
% Establish forecast error index
for j = 1:V.nfore
    eval(['V.fe' V.foretypes{j} ' = j;'])
end
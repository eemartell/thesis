function V = variables

% Variables: Names variables and assigns locations
% Inputs:
%     [None]
% Output:
%     V     :   Structure of variable locations and names

% Variables names
V.names = { 'n'          %1  Labor hours
            'w'          %2  Real wage rate
            'i'          %3  Nominal interest rate
            'in'         %4  Notional interest rate
            'c'          %5  Consumption
            'g'          %6  Growth rate
            's'          %7  Risk premium shock
            'y'          %8  Output
            'pi'         %9  Gross inflation
            'lam'        %11 Inverse MUC
            'mp'         %12 Monetary policy shock
          };

% Variables titles        
V.desc = {  'Labor Hours'
            'Real Wage Rate'
            'Nominal Interest Rate'  
            'Notional Interest Rate'
            'Consumption'
            'Growth'
            'Risk premium shock'
            'Output'           
            'Inflation Rate'
            'Inverse MUC'
            'Monetary Policy Shock'
         };

% Variables names for nonlinear model
V.plotnames = V.names;
V.nplotvar = length(V.plotnames);

% Variables descriptions for plotting
V.desc = V.desc;
      
% Shocks
V.shocktypes = {'g','s','mp'};  
V.nshock = length(V.shocktypes);

% Forecast errors
V.foretypes = {'lam','pi','g'}; 
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

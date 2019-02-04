function V = variables

% Variables: Names variables and assigns locations
% Inputs:
%     [None]
% Output:
%     V     :   Structure of variable locations and names

% Variables names
V.names = { 'c'          %1  Consumption            
            'i'          %2  Nominal interest rate
            'in'         %3  Notional interest rate
            'lam'        %4 Inverse MUC
            'w'          %5  Real wage rate
            'pi'         %6  Gross inflation
            'y'          %7  Output
            'n'          %8  Labor hours
            'g'          %9  Growth rate
            'a'          %10  Risk premium shock
            'mp'         %11 Monetary policy shock            
          };

% Variables titles        
V.desc = {  'Consumption'
            'Nominal Interest Rate'  
            'Notional Interest Rate'
            'Inverse MUC'
            'Real Wage Rate'
            'Inflation Rate'
            'Output'
            'Labor Hours'
            'Growth'
            'Preference shock'
            'Monetary Policy Shock'
         };

% Variables names for nonlinear model
V.plotnames = V.names;
V.nplotvar = length(V.plotnames);

% Variables descriptions for plotting
V.desc = V.desc;
      
% Shocks
V.shocktypes = {'g','a','mp'};  
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

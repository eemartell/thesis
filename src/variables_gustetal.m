function V = variables_gustetal

% Variables: Names variables and assigns locations
% Inputs:
%     [None]
% Output:
%     V     :   Structure of variable locations and names

% Variables names
V.names = { 'c'          %1  Consumption, non-ZLB      
            'c_zlb'      %2  Consumption, ZLB
            'i'          %3  Nominal interest rate
            'in'         %4  Notional interest rate
            'lam'        %5  Inverse MUC
            'w'          %6  Real wage rate
            'pi'         %7  Gross inflation
            'y'          %8  Output
            'n'          %9  Labor hours
            'g'          %10  Growth rate
            's'          %11  Risk premium shock
            'mp'         %12 Monetary policy shock            
          };

% Variables titles        
V.desc = {  'Consumption, non-ZLB'
            'Consumption, ZLB'
            'Nominal Interest Rate'  
            'Notional Interest Rate'
            'Inverse MUC'
            'Real Wage Rate'
            'Inflation Rate'
            'Output'
            'Labor Hours'
            'Growth'
            'Risk premium shock'
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

function V = variables

% Variables: Names variables and assigns locations
% Inputs:
%     [None]
% Output:
%     V     :   Structure of variable locations and names

% Variables names
V.names = { 'n'          %1  Labor hours
            'w'          %2  Real wage rate
            'c'          %3  Consumption
            'r'          %4  Nominal interest rate 
            'pi'         %5  Gross inflation
            'y'          %6  Output
            'psi'        %7  Marginal cost
            'z'          %8  Productivity
            'beta'       %9  Discount factor
            'k'          %10  Capital
            'i'          %11 Investment
            'rk'         %12 Rental Rate
            'q'          %13 Tobin's q
          };

% Variables names for plotting
V.plotnames = [V.names
               'realr'
               'rotcost'
               'ytil'
               'Ec'
               'Er'
               'Erk'
               'Ecper'];       
      
% Variables titles        
V.desc = {  'Labor Hours'
            'Real Wage Rate'
            'Consumption'
            'Nominal Interest Rate'
            'Inflation Rate'
            'Output'
            'Real Marginal Cost'
            'Technology'
            'Discount Factor'   
            'Capital'
            'Investment'
            'Real Rental Rate'
            'Tobin''s q'
            'Real Interest Rate'
            'Rot. Adjustment Cost'
            'Adjusted Output'
            'Expected Consumption'
            'Exp. Nom. Interest Rate'
            'Exp. Real Rental Rate'
            'Exp Cons. Growth'
         };

%Shocktypes: 
V.shocktypes = {'z','beta'};  
% Forecast errors
V.foretypes = {'c','pi','rk','i','q','beta'}; 
             
% Number of variables
V.nvar = length(V.names);
V.nplotvar = length(V.plotnames);
% Number of shocks
V.nshock = length(V.shocktypes);
% Number of forecast errors
V.nfore = length(V.foretypes);

% Establish variable index
for j = 1:V.nplotvar
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
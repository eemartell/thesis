function V = variables

% Variables: Names variables and assigns locations
% Inputs:
%     model
% Output:
%     V     :   Structure of variable locations and names

% Variables names
V.names =      {'n'          %1  Labor hours
                'w'          %2  Real wage rate
                'c'          %3  Consumption
                'y'          %4  Output                
                'k'          %5  Capital
                'i'          %6  Investment
                'q'          %7  Tobin's q 
                'rk'         %8  Rental rate
                'psi'        %9  Marginal cost 
                'r'          %10 Nominal interest rate
                'pi'         %11 Gross inflation
                'rn'         %12 Notional interest rate 
                'beta'       %13 Discount factor
                'g'          %14 Growth rate            
                'rgdp'       %15 Real GDP
                'clag'       %16 Lagged consumption
                'ilag'       %17 Lagged investment                
               };

% Variables titles        
V.desc =      { 'Labor Hours'
                'Real Wage Rate'
                'Consumption'
                'Output'
                'Capital'
                'Investment'
                'Tobin''s q'
                'Rental Rate'
                'Marginal Cost'
                'Nominal Interest Rate'
                'Inflation Rate'
                'Notional Interest Rate'
                'Discount Factor'               
                'Growth Rate'
                'Real GDP'
                'Lagged Consumption'
                'Lagged Investment'            
            };

%Shocktypes: 
V.shocktypes = {'g','beta','mp'};  
% Forecast errors
V.foretypes = {'c','beta','g','rk','q','i','pi'}; 

% Number of variables
V.nvar = length(V.names);
% Number of shocks
V.nshock = length(V.shocktypes);
% Number of forecast errors
V.nfore = length(V.foretypes);

% Establish variable index
for j = 1:V.nvar
   eval(['V.' V.names{j} ' = j;']);      
end
% Establish shock index
for j = 1:V.nshock
    eval(['V.eps' V.shocktypes{j} ' = j;'])
end
% Establish forecast error index
for j = 1:V.nfore
    eval(['V.fe' V.foretypes{j} ' = j;'])
end
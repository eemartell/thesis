function V = variables

% Variables: Names variables and assigns locations
% Inputs:
%     model
% Output:
%     V     :   Structure of variable locations and names

% Variables names
V.names =      {'c'          %1  Consumption 
                'n'          %2  Labor
                'x'          %3  Investment
                'k'          %4  Capital
                'y'          %5  Output
                'yf'         %7  Final goods       
                'u'          %5  Utilization cost
                'ups'        %8  Utilization choice
                'wg'         %9  Real wage growth gap    
                'xg'         %10 Investment growth
                'yg'         %11 Output growth
                'w'          %12 Real wage rate
                'wf'         %13 Flexible real wage
                'rk'         %14 Real rental rate
                'pi'         %15 Inflation rate
                'i'          %16 Nominal interest rate
                'in'         %17 Notional interest rate
                'q'          %18 Tobin's q
                'mc'         %19 Real marginal cost
                'lam'        %20 Inverse marginal utility of wealth  
                'g'          %21 Growth shock       
                's'          %22 Risk premium shock
                'mp'         %23 Monetary policy shock    
               };

% Variables titles        
V.desc =      { 'Consumption'
                'Labor'
                'Investment'
                'Capital'
                'Output'
                'Final Goods'
                'Utilization cost'
                'Utilization choice'
                'Real Wage Growth Gap'
                'Investment Growth Gap'
                'Output Growth Gap'
                'Real Wage'
                'Flexible Real Wage'
                'Real Rental Rate'
                'Inflation Rate'
                'Nominal Interest Rate'
                'Notional Interest Rate'
                'Tobin''s q'
                'Real Marginal Cost'
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
V.foretypes = {'lam','g','pi','rk','ups','q','xg','u','wg'}; 
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
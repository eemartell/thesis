clear all
close all
clc
%simulates 10,000 periods and stores Euler equation errors
%creates histogram of simulation and Euler equation errors

load('solutions/solutionti1.mat')

n = 10000;
e = normrnd(0,1,[n,1]);
u = normrnd(0,1,[n,1]);
v = normrnd(0,1,[n,1]);

sims = simulation_test(pf,P,S,G,V,e,u,v);

plotvars = [V.c, V.pi, V.i];
%   subplot padding
subpadbot = .295; % Increase if xaxis labels are cut off
subpadtop = .175; % Increase if subplot title is cut off
subpadleft = .1; % Increase if ylabel is cut off
subpadright = .05; % Decrease if white-space on RHS 

ylabtxt{1} = '\textbf{Consumption ($c$)}';
ylabtxt{2} = '\textbf{Inflation ($\pi$)}';
ylabtxt{3} = '\textbf{Interest rate ($i$)}';

legvar{1} = '$c$';
legvar{2} = '$\pi$';
legvar{3} = '$i$';

loc{1} = 'northwest';
loc{2} = 'northwest';
loc{3} = 'southeast';

figure('Color', 'white'); hold on; box on
plotdim = [length(plotvars) 1];

for i = 1:length(plotvars)
[col,row] = ind2sub(plotdim([2 1]),i);
left = (col-1+subpadleft)/plotdim(2);
bottom = 1-(row-subpadbot)/plotdim(1);
width = (1-(subpadleft+subpadright))/plotdim(2);
height = (1-subpadbot-subpadtop)/plotdim(1);
subplot(3,1,i, 'Position', [left bottom width height]); hold on; box on
histogram(sims(:,plotvars(i),:))
title(ylabtxt{i}, 'Interpreter', 'latex',...
      'FontName', 'Times New Roman')

legtext = ['10,000 realizations of ' ,legvar{i}];
set(legend(legtext),'Interpreter','latex',...
    'Location', loc{i}); 
  
end

resid{1} = log10(abs(sims(2:end,V.nplotvar+1,:)));
resid{2} = log10(abs(sims(2:end,V.nplotvar+2,:)));

ylabtxt{1} = 'Euler Equation'; 
ylabtxt{2} = 'Phillips Curve';

figure('Color', 'white'); hold on; box on
plotdim = [2 1];

for j = 1:2
[col,row] = ind2sub(plotdim([2 1]),j);
left = (col-1+subpadleft)/plotdim(2);
bottom = 1-(row-subpadbot)/plotdim(1);
width = (1-(subpadleft+subpadright))/plotdim(2);
height = (1-subpadbot-subpadtop)/plotdim(1);  
subplot(2,1,j, 'Position', [left bottom width height]); hold on; box on
histogram(resid{j})
title(ylabtxt{j}, 'Interpreter', 'latex',...
      'FontName', 'Times New Roman')
xlabel('Errors ($\log_{10}$)', 'Interpreter', 'latex')
end


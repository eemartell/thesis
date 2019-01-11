% Probabilities of going to the ZLB, given growth at 
% its steady state and the monetary policy shock at 
% its steady state 

% This plots the optimal choice of consumption as a 
% function of the demand shock

%Results: the plot goes below 1 in its last spot 

%TO DO: shade graph to see where it's binding 
%make my graph look like Throckmorton's plot

%%%calculate the nominal interest rate policy
r = max(1, inp);

%%%pick a cross section of the state space. look at the risk premium and
%%%set the other shock processes to steady state
aCS = find(r(ceil(O.g_pts/2),:,ceil(O.mp_pts/2))==1,1,'first');

%%%shade the ZLB region
aper = G.s_grid;%100*(G.a_grid-1);
firstZLB = find(r(:,aCS,:) == 1,1,'first');
shadeint = [aper(firstZLB) aper(end)];

x = G.s_grid;
y = pf.hh(ceil(O.g_pts/2),:,ceil(O.mp_pts/2)); %%only thing
figure('Color', 'white',...
       'pos', [8 8 800 350])
area(shadeint,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
    'HandleVisibility','off')
hold on
plot(x, y,... 
    'Color','k',...
    'Linewidth', 1.8)
xlabel('Risk Premium',...
        'FontName','Times New Roman')
ylabel('Consumption',...
        'FontName', 'Times New Roman')
%plot([G.s_grid(1) G.s_grid(11)] ,[S.c S.c], 'Color','k', 'Linewidth', 1)
xtks = 0.97:0.01:1.03;
xticks(xtks)
xlm = [x(1) x(end)];
xlim(xlm)
ytks = 0.328:0.001:0.336;
yticks(ytks)
ylm = [y(end) y(1)];
ylim(ylm)
for x = xtks
    plot([x x], ylm,...
        'LineStyle', ':',...
        'Color', 'k')
end
for y = ytks
    plot(xlim, [y y],...
        'LineStyle', ':',...
        'Color', 'k')
end
text(.993, .3338, '$\bar{\pi}$ = 0.3333', 'Interpreter', 'Latex', 'Fontsize', 10)
text(.993, .33355, '\downarrow', 'Fontsize', 10)

legtext = ['(\rho_a, \sigma_e, \sigma_u, \sigma_v) = (' ,num2str(P.rhoa),', ', num2str(P.sige),', ', num2str(P.sigu), ', ' num2str(P.sigv), ')'];
set(legend(legtext),'Interpreter','tex',...
    'Location', 'southeast'); 
hold off
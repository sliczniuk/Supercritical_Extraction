hold on
DATA_estimated = readtable('estimation.csv');
for ii=2:5
    KOUT = table2array(DATA_estimated(:,ii));
    plot(0:0.5:C0solid, KOUT(2) .* 10^-14 .* exp(KOUT(4) .* (1-(C0solid:-0.5:0)./C0solid) ) .*10^(12) ,'LineWidth', 4 )
end
hold off
legend({'$40[^\circ C] | 200 [bar] | 845[kg/m^3]$', '$50[^\circ C] | 200 [bar] | 795[kg/m^3]$', '$40[^\circ C] | 300 [bar] | 939[kg/m^3]$', '$50[^\circ C] | 300 [bar] | 903[kg/m^3]$'}, 'Location','northwest')
legend('boxoff')
xlabel('$c_s [kg/m^3]$')
ylabel('$D_i=(D_i^R \times \gamma) [m^2/s] \times 10^{-12}$')
fontsize(36,"points")
set(gcf,'PaperOrientation','landscape'); print(figure(1),['D_i_Uppsilon_results.pdf'],'-dpdf','-bestfit');
close all;
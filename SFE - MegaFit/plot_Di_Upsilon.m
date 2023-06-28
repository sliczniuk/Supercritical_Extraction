hold on
for ii=1:numel(DATA_set) 
    KOUT = DATA_K_OUT(:,ii);
    plot(0:0.5:C0solid, KOUT(1) .* 10^-14 .* exp(KOUT(4) .* (1-(C0solid:-0.5:0)./C0solid) ) ,'LineWidth',2 )
end
hold off
legend({'$40^\circ C / 200 bar$', '$50^\circ C / 200 bar$', '$40^\circ C / 300 bar$', '$50^\circ C / 300 bar$'}, 'Location','northwest')
legend('boxoff')
xlabel('$c_s [kg/m^3]$')
ylabel('$D_i \gamma [m^2/s]$')
fontsize(24,"points")
set(gcf,'PaperOrientation','landscape'); print(figure(1),['D_i_Uppsilon_results.pdf'],'-dpdf','-bestfit');
close all;
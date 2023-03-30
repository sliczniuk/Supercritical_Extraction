%RHO         = [840, 840, 785, 758, 910, 910, 870];
NAME        = {'$km[-]$', '$D_i[m^2/s]$', '$D_x[m^2/s]$', 'Toral mass of oil', 'Mass ratio', '$\sigma$'};
%EXP_Name    = {'exp1','exp2'};

%% Fit the parameters based as linear function of density

figure(1)

for jj=1:2
    for ii=1:Nk
        
        [trendLine, SLine]   = polyfit(RHO, DATA_K_OUT(ii,:), jj );
        p_est                = polyval(trendLine, RHO);
        p_y                  = polyval(trendLine, linspace(min(RHO), max(RHO) ) );
        
        subplot(2,3,ii)
        scatter(RHO, DATA_K_OUT(ii,:), 'filled')
        hold on
        plot(linspace(min(RHO), max(RHO) ),p_y,'LineWidth',2);
        hold off
    
        title(['$R^2 = $', mat2str( round(power(corr2(DATA_K_OUT(ii,:), p_est),2),2) ) ])
    
        xlabel('$\rho~[kg/m^3]$')
        ylabel(NAME{ii})
    
    end
    
    set(gcf,'PaperOrientation','landscape'); print(figure(1),['Trend_Lines_order_',mat2str(jj),'.pdf'],'-dpdf','-bestfit'); 
    close all;

end

%{
figure(1)
for jj=1:2
    for ii=1:3
        
        f = fit(RHO', DATA_K_OUT(ii,:)',EXP_Name{jj});
        
        subplot(3,1,ii)
        scatter(RHO, DATA_K_OUT(ii,:), 'filled')

        hold on
        plot(linspace(min(RHO), max(RHO) ),f(linspace(min(RHO), max(RHO) )),'LineWidth',2);
        hold off
    
        title(['$R^2 = $', mat2str( round(power(corr2(DATA_K_OUT(ii,:), f(RHO)' ),2),2) ) ])
    
        xlabel('$\rho~[kg/m^3]$')
        ylabel(NAME{ii})
    
    end
    
    set(gcf,'PaperOrientation','landscape'); print(figure(1),['Trend_EXP_',mat2str(jj),'.pdf'],'-dpdf','-bestfit'); close all;

end
%}
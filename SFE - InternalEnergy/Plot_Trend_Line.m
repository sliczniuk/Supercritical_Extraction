RHO         = [840, 785, 910, 870];
NAME        = {'$km[-]$', '$D_i[m^2/s]$', '$D_x[m^2/s]$', '$C_{sat}[kg/m^3]$', 'mSOL_s', 'max_ms'};

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
    
    set(gcf,'PaperOrientation','landscape'); print(figure(1),['Inital_state_2_Trend_Lines_order_',mat2str(jj),'_No_Delay.pdf'],'-dpdf','-bestfit'); close all;

end
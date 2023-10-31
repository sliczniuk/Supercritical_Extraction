%RHO         = [840, 840, 785, 758, 910, 910, 870];
%NAME        = {'$k_m[-]$', '$D_i[m^2/s]$', '$D_e^M[m^2/s]$', '$m_{total}$', '$\tau$', '$\sigma$'};
%EXP_Name    = {'exp1','exp2'};

%% Fit the parameters based as linear function of density

figure(1)

for jj=1:2
    for ii=1:Nk
        
        [trendLine, SLine]   = polyfit(RHO, DATA_K_OUT(ii,:), jj );
        p_est                = polyval(trendLine, RHO);
        p_y                  = polyval(trendLine, linspace(min(RHO), max(RHO) ) );
        
        subplot(3,3,ii)
        scatter(RHO, DATA_K_OUT(ii,:), 'filled')
        hold on
        plot(linspace(min(RHO), max(RHO) ),p_y,'LineWidth',2);
        hold off

        if jj==1
            pp=sprintf('%s=$(%6.4g) \\rho + (%6.4g)$\n',NAME{ii},trendLine');
        else
            pp=sprintf('%s=$(%6.4g) \\rho^2 + (%6.4g) \\rho + (%6.4g)$\n',NAME{ii},trendLine');
        end

        oo = sprintf('$R^2 = %g$\n', round(power(corr2(DATA_K_OUT(ii,:), p_est),2),2) ) ;
    
        title( [oo,pp], 'FontSize', 8)
    
        xlabel('$\rho~[kg/m^3]$','Interpreter','latex')
        ylabel(NAME{ii},'Interpreter','latex')
    
    end
    
    set(gcf,'PaperOrientation','landscape'); print(figure(1),['Trend_Lines_order_',mat2str(jj),'.pdf'],'-dpdf','-bestfit'); 
    close all;

end

for jj=1:2
    for ii=1:Nk
        
        [trendLine, SLine]   = polyfit(RE, DATA_K_OUT(ii,:), jj );
        p_est                = polyval(trendLine, RE);
        p_y                  = polyval(trendLine, linspace(min(RE), max(RE) ) );
        
        subplot(3,3,ii)
        scatter(RE, DATA_K_OUT(ii,:), 'filled')
        hold on
        plot(linspace(min(RE), max(RE) ),p_y,'LineWidth',2);
        hold off

        if jj==1
            pp=sprintf('%s=$(%6.4g) Re + (%6.4g)$\n',NAME{ii},trendLine');
        else
            pp=sprintf('%s=$(%6.4g) Re^2 + (%6.4g) Re + (%6.4g)$\n',NAME{ii},trendLine');
        end

        oo = sprintf('$R^2 = %g$\n', round(power(corr2(DATA_K_OUT(ii,:), p_est),2),2) ) ;
    
        title( [oo,pp], 'FontSize', 8)
    
        xlabel('$Reynolds~number$','Interpreter','latex')
        ylabel(NAME{ii},'Interpreter','latex')
    
    end
    
    set(gcf,'PaperOrientation','landscape'); print(figure(1),['Trend_Lines_RE_order_',mat2str(jj),'.pdf'],'-dpdf','-bestfit'); 
    close all;

end
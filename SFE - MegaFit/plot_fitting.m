function [] = plot_fitting(DATA, Time_in_sec, RHO, F, g, x0, KOUT, which_parameter, N_Delay, N_Sample, theta, print_mode)

    Time = 0:numel(Time_in_sec);
    
    for i = 1:numel(DATA)
    
        subplot(2,2,i)

        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)], ...
            'renderer','OpenGL');
        set(gca,'TickLabelInterpreter','latex')

        load(DATA{i});

        %T0homog   = 50 + 273.15;                             % Extractor initial temperature (pseudo-homogeneous)
        feedTemp  = T0homog      * ones(1,length(Time_in_sec));  % Kelvin
        feedPress = feedPress(1) * ones(1,length(Time_in_sec));  % Bars

        feedFlow  = 0.4 * RHO(i) * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec
        
        %% Inital values

        uu = [feedTemp', feedPress', feedFlow'];

        for k=1:1%numel(DATA)
    
            %[uu(1,:), KOUT(:,k)']
    
            %% Simulate system with estiamted parameter
            [yy_out, tt_out, xx_out] = simulateSystem_ParameterEstimation(F, g, x0, uu, KOUT(:,i));
    
            yy_out_D = [ zeros(1,N_Delay) yy_out(1:end-N_Delay)];

            %legend_name = sprintf('km = %g, Di = %g, \\sigma = %g', KOUT(:,k)');

            hold on
%            plot(Time, yy_out_D, 'LineWidth',1, 'DisplayName',legend_name);
            %plot(Time(1:N_Sample:end), data, 'o','MarkerSize',2)
            plot(Time, yy_out_D, 'LineWidth',1);
            plot(Time(1:N_Sample:end), data ,'o')
            hold off
            
        end

        xlim([0,Time(end)])
        ylim([-10,110])

        yline(50,'--')
        xline(50,'--')

        yline(75,'--')
        xline(100,'--')

        %caption = sprintf('T = %g [K], P = %g [bar]\n km = %g, Di = %g, Dx = %g, $\\sigma$ = %g', T0homog, feedPress(1), KOUT(:,i)');
        caption = sprintf('T = %g [K], P = %g [bar]', feedTemp(1), feedPress(1));
        %title(caption, 'FontSize', 4, 'interpreter','latex');
        title(caption, 'interpreter','latex');

        xlabel('Time [min]','interpreter','latex')
        ylabel('Yield [\%]','interpreter','latex')
    
    end
    %% save data
    if (print_mode=='print_on')
        print(gcf, 'fitting.pdf','-dpdf','-r700','-bestfit'); close all
    end

end
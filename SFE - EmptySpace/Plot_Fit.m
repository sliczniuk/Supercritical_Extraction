    %figure('Visible','off')
    figure(ii)

    h = tiledlayout(2,2);

    nexttile
    %pbaspect([3 1 1])
    hold on
%    plot(Time, xx_0(end,:));
    plot(Time, xx_out(end,:));
    plot(SAMPLE, data_org,'o');
    xline(PreparationTime,'--k')
    ylim([0 80])
    hold off
    legend('$k_0$', '$k_{out}$', 'data', 'Location','northoutside')
    legend('Orientation','horizontal')
    legend('boxoff')
    ylabel('CFD: y(t)','Interpreter','latex')
    xlabel('Time [min]','Interpreter','latex')

    nexttile
    %pbaspect([3 1 1])
    hold on
%    plot(SAMPLE(1:end-1), diff(xx_0(end,N_Sample))   );
    plot(SAMPLE(1:end-1), diff(xx_out(end,N_Sample)) );
    plot(SAMPLE(1:end-1), diff(data_org),'o');
    xline(PreparationTime,'--k')
    hold off
    ylim([0 15])
    legend('$k_0$', '$k_{out}$', 'data', 'Location','northoutside')
    legend('Orientation','horizontal')
    legend('boxoff')
    ylabel('PDF: norm $\frac{dy(t)}{dt}$','Interpreter','latex')
    xlabel('Time [min]','Interpreter','latex')

    nexttile
    %pbaspect([3 1 1])
    hold on
    %plot(Time, 1e3 * (sum(xx_0(  0*nstages+1:1*nstages,:) .* V_fluid) + sum(xx_0(  1*nstages+1:2*nstages,:) .* V_solid)) + xx_0(  Nx,:))
    plot(Time, 1e3 * (sum(xx_out(0*nstages+1:1*nstages,:) .* V_fluid) + sum(xx_out(1*nstages+1:2*nstages,:) .* V_solid)) + xx_out(Nx,:))
    hold off
    ylabel('Total mass of solute [g]','Interpreter','latex')
    xlabel('Time [min]','Interpreter','latex')

    nexttile
    %pbaspect([3 1 1])
    hold on
    %plot(Time, 1e3 * (sum(xx_0(  0*nstages+1:1*nstages,:) .* V_fluid) + sum(xx_0(  1*nstages+1:2*nstages,:) .* V_solid)) + xx_0(  Nx,:))
    plot(Time, [feedFlow, feedFlow(end)] )
    hold off
    ylabel('Superficial velocity [m/s]','Interpreter','latex')
    xlabel('Time [min]','Interpreter','latex')

    subtitle([DATA,'_',mat2str(round([DATA_K_OUT(:,ii)],3))],'Interpreter','latex')

    set(gcf,'PaperOrientation','landscape'); print(figure(ii),['Yield_',DATA,'_F',mat2str(V_Flow),'_No_Delay.pdf'],'-dpdf','-bestfit'); close all

    %%
%{
    figure('Visible','off')
    NAME = {'c_f','c_s','T','\rho'};

    h = tiledlayout(2,2);

    for i=1:4

        nexttile
        imagesc(Time,L_nstages./L,xx_out((i-1)*nstages+1:i*nstages,:)); colormap jet; colorbar
        hold on
        yline(L_nstages(nstagesbed([1,end]))./L,'--w')
        xline(PreparationTime,'--w')
        hold off
        pbaspect([2 1 1])
        c = colorbar;
        title(['$',NAME{i},'$'],'Interpreter','latex')
        set(c,'TickLabelInterpreter','latex')
        ylabel('$L [-]$','Interpreter','latex')
        xlabel('Time [min]','Interpreter','latex')
        %axis square
    end

    set(gcf,'PaperOrientation','landscape'); print(figure,['Profiles_',DATA,'_F',mat2str(V_Flow),'.pdf'],'-dpdf','-bestfit'); close all
%}
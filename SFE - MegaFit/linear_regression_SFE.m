function [Results] = linear_regression_SFE(KOUT,RHO,Operating_Conditions_experiments, print_mode)
    
    dataTemp   = Operating_Conditions_experiments(1,:);
    dataPress  = Operating_Conditions_experiments(2,:);
    
    [X,Y] = meshgrid(min(dataTemp):max(dataTemp),min(dataPress):max(dataPress));
    
    figure()
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)], ...
        'renderer','OpenGL');
    
    Name = {'km','Di','Dx'};

    Results = {};
    
    for i=1:size(KOUT,1)-1
    
        % linear regression - rho
        data = KOUT(i,:);
        [A,B] = fit(RHO',data','poly1');
        subplot(2,size(KOUT,1)-1,i)
        set(gca,'TickLabelInterpreter','latex')
        hold on
        plot(RHO, A.p1*RHO+A.p2)
        plot(RHO, data, '.')
        ylabel(Name{i},'interpreter','latex')
        xlabel('RHO [kg/m3]','interpreter','latex')
        str = sprintf('$R^2 =$%.3f',B.rsquare);
        title(str,'interpreter','latex')
        hold off

        % linear regression - T,P plane
        [AA,BB] = fit([dataTemp',dataPress'],data','poly11');
        Z = AA.p00 + AA.p10* X +AA.p01 * Y;
        subplot(2,size(KOUT,1)-1,size(KOUT,1)-1+i)
        set(gca,'TickLabelInterpreter','latex')
        hold on
        surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.5)
        plot3(dataTemp,dataPress,data,'.')
        hold off
        ylabel('P [bar] ','interpreter','latex')
        xlabel('T [K]','interpreter','latex')
        zlabel(Name{i},'interpreter','latex')
        str = sprintf('$R^2 =$%.3f',BB.rsquare);
        title(str,'interpreter','latex')
        view(45,45)

        Results{i} = AA;

    end
    
    %% save data
    if (print_mode=='print_on')
        print(gcf, 'fitting_regression.pdf','-dpdf','-r700'); close all
    end

end
function [EA, km0] = Arrhenius_km(T,kp,R,rho_s)
    
    rho = [485.15 285 235.4	212.4];
    
    km = kp.*rho / rho_s;

    [P,S] = polyfit(1./T,log(km),1);

    EA = -P(1)*R;
    km0 = exp(P(2));
    
%     TT = linspace(T(1),T(end));
%     hold on
%     plot(1./T,log(km),'o')
%     plot(1./TT, P(1).*(1./TT) + P(2),'-')
%     xlabel('1/T [1/K]','interpreter','latex')
%     ylabel('$ln(km) [-]$','interpreter','latex')
%     R2 = 1 - (S.normr/norm(log(km) - mean(log(km))))^2;
%     str = sprintf('$R^2 =$%.3f',R2);
%     annotation('textbox',[.8 .6 .3 .3],'string',str,'FitBoxToText','on','EdgeColor','black','FontSize',15,'interpreter','latex')   
%     hold off
%     title(['Line Equation is y = F(x) = ',num2str(P(1)),'*x + ','(',num2str(P(2)),')'],'interpreter','latex');
    
end
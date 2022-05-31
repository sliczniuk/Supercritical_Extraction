function [EA, D0] = Arrhenius_Di(T,D, R)

    [P,S] = polyfit(1./T,log(D),1);

    EA = -P(1)*R;
    D0 = exp(P(2));
    
    %{
    dT  = [ ones(numel(T),1) 1./T' ];        % Data temperatures of some biomass
    EA1 = -R * log10(D(1) / D(end)) / (dT(1,2) - dT(numel(T),2)); % Activation energies
                                                        % We assume that dDi = beta * dT and estimate beta using least-squares
    betah1 = inv(dT'*dT)*dT'*D';                            % Estimated coefficients
    betah = betah1(1);                                       %

    
    %EA1 = -R * log(D(1) / D(end)) / (dT(1,2) - dT(4,2)); % Activation energies
                                                        % We assume that dDi = beta * dT and estimate beta using least-squares
    %betah = inv(dT'*dT)*dT'*D';                            % Estimated coefficients
    %betah = betah(1); 
    %}
%     TT = linspace(T(1),T(end));
%     hold on
%     plot(1./T,log(D),'o')
%     plot(1./TT, P(1).*(1./TT) + P(2),'-')
%     xlabel('1/T [1/K]','interpreter','latex')
%     ylabel('$ln(Di) [m^2/s]$','interpreter','latex')
%     R2 = 1 - (S.normr/norm(log(D) - mean(log(D))))^2;
%     str = sprintf('$R^2 =$%.3f',R2);
%     annotation('textbox',[.8 .6 .3 .3],'string',str,'FitBoxToText','on','EdgeColor','black','FontSize',15,'interpreter','latex')   
%     hold off
%     title(['Line Equation is y = F(x) = ',num2str(P(1)),'*x + ','(',num2str(P(2)),')'],'interpreter','latex');
%     
end
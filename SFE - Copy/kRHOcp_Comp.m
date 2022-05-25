function krhocp = kRHOcp_Comp(T,P,Z,RHO,CP,theta)

    % http://www.ffrc.fi/FlameDays_2009/3B/HankalinPaper.pdf
    
    epsi     = theta{4};     % Void bed fraction    
    rhoSolid = theta{7};     % Density of solid
    cpSolid  = theta{24};      % 1.5* 10^3;  kJ/K/ Kg -> J / K / Kg 
    
    k = HeatConductivity_Comp(T,P,Z,RHO,theta)   * 10^-3;        % mili W/m/K -> W/m/K
    cpFluid = CP;                % J/kg/K
    
    krhocp = k ./ (cpFluid.*(1-epsi).*RHO + cpSolid.*epsi*rhoSolid);

end

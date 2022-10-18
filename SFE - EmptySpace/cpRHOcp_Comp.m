function cprhocp = cpRHOcp_Comp(T, P, Z, RHO, CP, epsi, theta)

    %epsi     = theta{4};     % Void bed fraction    
    rhoSolid = theta{7};    % Density of solid
    cpSolid  = theta{24};   % [kJ/kg/K]
    
    cpFluid = CP;           % [kJ/kg/K]
    
    %cpSolid = 1.5E3; %theta(); % J / K / Kg 

   cprhocp = cpFluid ./ (cpFluid .* (1-epsi) .* RHO + cpSolid .* epsi .* rhoSolid);

end

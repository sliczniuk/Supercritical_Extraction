function drhodp_T = drhodp(T, P, RHO, theta)

    %% Parameters
    Tc      = theta{10};         % Critical temperature [K]
    Pc      = theta{11};         % Critical pressure [bar]
    R       = theta{12};         % Universal gas constant, [cm3-bar/K-mol]
    kappa   = theta{13};
    MW      = theta{14};         % Molar mass [g/mol]       

    %% Recalculate the Peng Robinson parameters

    RHO = RHO*1e-3;       %  kg/m3 -> g/cm3
    
    Tr = T ./ Tc;
    alpha = (1 + kappa .* (1 - sqrt(Tr))).^2;

    a = 0.4572350 .* R.^2 .* Tc.^2 ./ Pc;
    b = 0.0777961 .* R    .* Tc    ./ Pc;

    A = a .* alpha .* P ./ R.^2 ./ T.^2;
    B = b .* P ./ R ./ T;

    v = 1./RHO*MW;

    %%

    dPdV = -R.*T./(v-b).^2 + 2.*a.*(v+b)./(v.*(v+b) + b.*(v-b)).^2;

    dVdP = 1/dPdV ; 

    drhodp_T =  -1/(v^2) * dVdP * RHO * 1e3;
end
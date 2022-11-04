function beta = compressibility_beta(T,P,RHO,theta)
    
    %% Parameters
    Tc      = theta{10};         % Critical temperature [K]
    Pc      = theta{11};         % Critical pressure [bar]
    R       = 83.1447;    % Universal gas constant, [m3-bar/K-mol]
    kappa   = 0.2250;
    MW      = theta{14};         % Molar mass [g/mol]
    CP_0    = theta{17};
    CP_A    = theta{18};
    CP_B    = theta{19};
    CP_C    = theta{20};
    CP_D    = theta{21};
    
    %% Recalculate the Peng Robinson parameters
    Tr = T ./ Tc;
    Pr = P ./ Pc;
    %alpha = (1 + kappa .* (1 - sqrt(Tr))).^2;
    
    m     = 0.37464 + 1.54226*kappa - 0.26992 * kappa^2;
    alpha = ( 1 + m .* ( 1 - Tr.^0.5 ) ).^2;
    
    a = 0.45724 .* R.^2 .* Tc.^2 *alpha ./ Pc;
    b = 0.0777961 .* R    .* Tc           ./ Pc;
    
    %A = a .* alpha .* P ./ R.^2 ./ T.^2;
    A = 0.45723553.*alpha.*Pr./Tr.^2;
    B = b .* P ./ R ./ T;
    
    v = 1./RHO*MW*1e6;
    
    %% Derivatives
    
    %  (δP/ δV)_T = -RT/(v - b)² + 2a(v + b)/[v(v + b) + b(v - b)]²
    dPdV = -R.*T./(v-b).^2 + 2.*a.*(v+b)./(v.*(v+b) + b.*(v-b)).^2;
    

    %%

    dVdP = 1/dPdV;

    beta = -1./v .* dVdP ;

end
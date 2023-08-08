function alpha = ThermalExpansion(T, P, Z, RHO,  theta)

    %% Parameters
    Tc      = theta{10};         % Critical temperature [K]
    Pc      = theta{11};         % Critical pressure [bar]
    R       = 83.1447;    % Universal gas constant, [m3-bar/K-mol]
    kappa   = 0.2250;
    MW      = theta{14};         % Molar mass [g/mol]       
    
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
    
    %  (δa/ δT)_V = -kappa/[ ((T*Tc)^0.5) (1 + kappa( 1 - (T/Tc)^0.5))]
    %dadT = -kappa.*a./( sqrt(T*Tc).*(1 + kappa.*( 1 - sqrt(Tr))) );
    dadT = -m.*a./( sqrt(T*Tc).*(1 + m.*( 1 - sqrt(Tr))) );
    
    %  (δA/δT)_P = (P/(RT)²)(a' - 2a/T)
    dAdT = (P./(R.*T).^2).*(dadT - 2*a./T);

    %  (δB/δT)_P = -bP/(RT²)
    dBdT = -b.*P./(R.*T.^2);
    
    %  (δZ/ δT)_P = Num / Denom
    %   Num = (δA/δT)P (B-Z) + (δB/δT)_P (6BZ+2Z-3B²-2B+A-Z²)
    %   Denom = 3Z² + 2(B-1)Z + (A-2B-3B²)
    Num = dAdT.*(B-Z) + dBdT.*(6*B.*Z+2.*Z-3*B.^2-2*B+A-Z.^2);
    Denom = 3.*Z.^2 + 2.*(B-1).*Z + (A-2.*B-3.*B.^2);
    dZdT = Num./Denom;

    alpha = R./(P.*v) .* ( Z + T .* dZdT );

end
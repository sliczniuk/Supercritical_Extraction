function Cp = SpecificHeatComp(T, P, Z, RHO, theta)

    % Specific heat equations
    % https://cheguide.com/specific_heat_ratio.html
   
    % The article with specific heat plots, to compare
    % https://www.sciencedirect.com/science/article/abs/pii/S0017931014006577

    %% Parameters
    Tc      = theta{10};         % Critical temperature [K]
    Pc      = theta{11};         % Critical pressure [bar]
    R       = theta{12};    % Universal gas constant, [m3-bar/K-mol]
    kappa   = theta{13};
    MW      = theta{14};         % Molar mass [g/mol]       
    CP_0    = theta{17};     
    CP_A    = theta{18};       
    CP_B    = theta{19};       
    CP_C    = theta{20};       
    CP_D    = theta{21};

    RHO = RHO*1e-3;       %  kg/m3 -> g/cm3
    
    %% Recalculate the Peng Robinson parameters
    Tr = T ./ Tc;
    alpha = (1 + kappa .* (1 - sqrt(Tr))).^2;

    a = 0.4572350 .* R.^2 .* Tc.^2 ./ Pc;
    b = 0.0777961 .* R    .* Tc    ./ Pc;

    A = a .* alpha .* P ./ R.^2 ./ T.^2;
    B = b .* P ./ R ./ T;

    v = 1./RHO*MW;
    
    %% Derivatives
    
    %  (δP/ δV)_T = -RT/(v - b)² + 2a(v + b)/[v(v + b) + b(v - b)]²
    dPdV = -R.*T./(v-b).^2 + 2.*a.*(v+b)./(v.*(v+b) + b.*(v-b)).^2;
    
    %  (δa/ δT)_V = -kappa/[ ((T*Tc)^0.5) (1 + kappa( 1 - (T/Tc)^0.5))]
    dadT = -kappa.*a./( sqrt(T*Tc).*(1 + kappa.*( 1 - sqrt(Tr))) );
    
    %  (δA/δT)_P = (P/(RT)²)(a' - 2a/T)
    dAdT = (P./(R.*T).^2).*(dadT - 2*a./T);

    %  (δB/δT)_P = -bP/(RT²)
    dBdT = -b.*P./(R.*T.^2);
    
    %  (δP/ δT)V = R/(v - b) - a'/[v(v + b) + b(v - b)]
    dPdT = R./(v-b) - dadT./( v.*(v+b) + b.*(v-b) );
    
    %  (δZ/ δT)_P = Num / Denom
    %   Num = (δA/δT)P (B-Z) + (δB/δT)_P (6BZ+2Z-3B²-2B+A-Z²)
    %   Denom = 3Z² + 2(B-1)Z + (A-2B-3B²)
    Num = dAdT.*(B-Z) + dBdT.*(6*B.*Z+2.*Z-3*B.^2-2*B+A-Z.^2);
    Denom = 3.*Z.^2 + 2.*(B-1).*Z + (A-2.*B-3.*B.^2);
    dZdT = Num./Denom;
    
    %  (δV/ δT)_P = (R/P)[ T(δZ/δT)_P + Z]
    dVdT = (R./P) .* (T.*dZdT + Z);
    
    %  a" = a kappa (1 + kappa)(Tc/T)^0.5 / (2*T*Tc)
    d2adT2 = a.*kappa.*(1+kappa).*sqrt(Tc./T)./(2.*T.*Tc);

    %% Heat Capacity - Ideal gas
    Cp_Ideal = CP_0 * (CP_A + CP_B.*T + CP_C.*T.^2 + CP_D.*T.^3);
    Cv_Ideal = Cp_Ideal - R/10;

    %% Heat Capacity - correction for real gasses
    % CvR = (δUR/δT)V
    % UR = [(Ta'-a)/b(8)^0.5] ln[(Z+B(1+2^0.5))/(Z+B(1-2^0.5))]
    % CvR = [Ta"/b(8)^0.5] ln[(Z+B(1+2^0.5))/(Z+B(1-2^0.5))]
    Cv_corr = ((T.*d2adT2)./(b.*sqrt(8))) .* log( ( Z+B.*(1+sqrt(2)) )./(Z+B.*(1-sqrt(2))) )./10;
    
    % CpR = CvR + T*(δP/δT)_V (δV/δT)_P - R
    Cp_corr = Cv_corr + T.*dPdT.*dVdT/10 - R/10;

    %% Heat Capacity - Real gas
    Cv_mol = Cv_Ideal + Cv_corr;               %[J/mol/K]
    Cp_mol = Cp_Ideal + Cp_corr;               %[J/mol/K]

    Cv = Cv_mol/MW;                            %[J/g/K] = [kJ/kg/K]
    Cp = Cp_mol/MW;                            %[J/g/K] = [kJ/kg/K]

end
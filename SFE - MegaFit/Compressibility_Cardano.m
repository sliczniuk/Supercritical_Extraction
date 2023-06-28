function Z = Compressibility_Cardano(T,P,theta)
    %% Parameters
    Tc      = theta{10};         % Critical temperature [K]
    Pc      = theta{11};         % Critical pressure [bar]
    R       = 83.1447;           % Universal gas constant, [m3-bar/K-mol]
    kappa   = 0.2250;
    MW      = theta{14};         % Molar mass [g/mol]       
    
    CPA =  1.98E+01;
    CPB =  7.34E-02;
    CPC = -5.60E-05;
    CPD =  1.72E-08;

    %% Reference state %TODO: Add a function to calulate the properties of the reference state
    TREF = 298.15;   % K
    PREF = 0.101325; % MPa

    Delta_CP_IG_REF = -41.81755442;
    
    %% Recalculate the Peng Robinson parameters
    Tr = T ./ Tc;
    Pr = P ./ Pc;
    %alpha = (1 + kappa .* (1 - sqrt(Tr))).^2;
    
    m     = 0.37464 + 1.54226*kappa - 0.26992 * kappa^2;
    alpha = ( 1 + m .* ( 1 - Tr.^0.5 ) ).^2;

    a = 0.45724 .* R.^2 .* Tc.^2 *alpha ./ Pc;
    b = 0.0777961 .* R    .* Tc         ./ Pc;

    %A = a .* alpha .* P ./ R.^2 ./ T.^2;
    A = 0.45723553.*alpha.*Pr./Tr.^2;
    B = b .* P ./ R ./ T;

    UC = - (1 - B)              ;
    SC =   (A - 2.*B - 3.*B.^2) ;
    TC = - ( A .* B - B.^2 - B.^3);

    PC = (3.*SC - UC.^2) ./ 3;
    QC = (2.*UC.^3 ./ 27) - UC.*SC./3 + TC;

    DC = (PC./3).^3 + (QC./2).^2;
    
end
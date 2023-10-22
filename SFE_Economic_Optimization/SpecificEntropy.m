function [s_kg] = SpecificEntropy(T, P, Z, RHO, theta)

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

    %% Reference state
    TREF = 298.15;   % K
    PREF = 0.101325; % MPa

    Delta_S_IG_REF = -0.094155634;
    
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

    v = 1./RHO*MW*1e6;
    
    %% Derivatives
    
    Delta_S_IG_R = CPA.*log(T./TREF) + CPB.*(T-TREF) + CPC./2*(T.^2-TREF.^2) + CPD./3*(T.^3-TREF.^3) - R./10.*log(P/10./PREF);

    Delta_S_IG   = (R/10).*log(Z-B)-A.*(R/10)./B./2.8284.*m.*sqrt(Tr./alpha) .* log( (Z+2.4142.*B)./(Z-0.4142.*B) );

    s_mol         =   ( Delta_S_IG + Delta_S_IG_R - Delta_S_IG_REF );          % J/mol % minus sign is used to work with positive values
    s_kg          =   s_mol./MW;                                                   % J/mol -> J/kg
    s_kg          =   s_kg .* 1e-3;                                                % J/kg -> KJ/kg

end
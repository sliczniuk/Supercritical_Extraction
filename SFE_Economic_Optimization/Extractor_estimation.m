function [Cost] = Extractor_estimation(P,Parameters)

    r       = Parameters{3} .* 39.37;                                       % Radius of the extractor  [m->inch]
    L       = Parameters{6} .* 39.37;                                       % Total length of the extractor [m->inch]
    %r            = 78/2;
    %L            = 480;
    Max_Stress   = 13100;
    Weld_E       = 1;
    material_rho = 0.284;                                                   % carbon steel: 490 lb∕ft3 or 0.284 lb∕in3.
    FM           = 1.2;                                                     % Low-alloy steel - 1.2

    P_psig  = (P-1) .* 14.5037                                             % bar -> psig
    
    Pd      = exp(0.60608 + 0.91615*(log(P_psig)) + 0.0015655*(log(P_psig)).^2);
    tp      = Pd.*(2.*r) / (2 * Max_Stress * Weld_E - 1.2*P_psig);
    ts      = tp;
    W       = pi*(2.*r + ts) * (L + 0.8*(2*r)) * ts * material_rho;
    Cv      = exp(5.6336 + 0.4599*log(W) + 0.00582*log(W).^2);
    C       = FM * Cv;                                                      % account for a material factor
    Cost    = (600/567) * C;
end
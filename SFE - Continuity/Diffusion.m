function D = Diffusion(T,P,theta)

    MW     = theta{14}*10^3;
    M      = theta{21};
    epsi_p = theta{22};
    V      = 339709.057834;                     % the solute molar volume at boiling point (m3/kmol).
    mu     = Viscosity(T,P)*10^(-3);
    
    D12  = M .* sqrt(MW).*T./mu./V.^(0.6);
    D    = D12.*epsi_p./(2-epsi_p);
    
end
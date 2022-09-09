function v = Velocity(F,rho,theta)

    V       = theta{3};     % Extractor volume (m3)
    epsi    = theta{4};     % Void bed fraction
    L       = theta{6};     % Length of the extractor (m)

    A = V/L;
    
    v = F./(epsi * A .* rho);

end
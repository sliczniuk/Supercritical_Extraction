function v = Velocity(F,rho,epsi,theta)

    r = theta{3};     % Extractor length (m)

    A = pi*r^2 ;       % Cross-section of the extractor (m^2)
    
    v = F./(epsi .* A .* rho);

end
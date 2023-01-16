function D_L = axial_diffusion(T,P,F,rho,theta)

    V       = theta{3};     % Extractor volume (m3)
    epsi    = theta{4};     % Void bed fraction
    dp      = theta{5};
    L       = theta{6};     % Length of the extractor (m)
    
    a = theta{25}; % 0.1
    b = theta{26}; % 0.011
    c = theta{27}; % 0.48

    v       = Velocity(F,rho,theta);
    tau     = V/(F/rho);
    mu      = Viscosity(T,P,theta)./rho;
    
    Re = dp * v .* rho/epsi./mu;
    
    Pe = (a/epsi) + (b/epsi)*(epsi*Re).^(c);
    
    D_L = dp*v./epsi.*Pe;

end
function D_L = axial_diffusion(T, epsi, v, rho)

    mu = Viscosity(T,rho);

    dp = 0.15;

    Re = dp .* v .* rho .* epsi ./ mu;
    
    D_L = 0.07334 .* Re - 14.84;

end
function D = Diffusion(rho, parameters)

    a = parameters{44};
    b = parameters{45};

    D =  a.* rho + b;
    
end
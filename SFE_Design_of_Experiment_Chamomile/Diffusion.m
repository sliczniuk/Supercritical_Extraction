function D = Diffusion(Re, F, parameters)

    a = parameters{44};
    b = parameters{45};

    D =  -0.131 - 8.337 .* Re + 0.687 * (F * 10^5);
    
end
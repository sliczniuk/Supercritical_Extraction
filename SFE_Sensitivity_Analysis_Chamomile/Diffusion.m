function D = Diffusion(RE, parameters)

    a = parameters{44};
    b = parameters{45};

    D =  a.* RE + b;
    
end
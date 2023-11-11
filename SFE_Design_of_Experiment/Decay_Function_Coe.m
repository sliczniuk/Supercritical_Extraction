function [Upsilon] = Decay_Function_Coe(rho, parameters)

    a = parameters{46};
    b = parameters{47};

    Upsilon  = a .* rho + b;

end
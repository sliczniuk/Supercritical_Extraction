function [Upsilon] = Decay_Function_Coe(Re, parameters)

    a = parameters{46};
    b = parameters{47};

    Upsilon  = a .* Re + b;

end
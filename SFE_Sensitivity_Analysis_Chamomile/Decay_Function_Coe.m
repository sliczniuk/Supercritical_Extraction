function [Upsilon] = Decay_Function_Coe(Re, F, parameters)

    a = parameters{46};
    b = parameters{47};

    Upsilon  = 4.717 + 11.012 .* Re - 0.876 .* (F * 10^5);

end
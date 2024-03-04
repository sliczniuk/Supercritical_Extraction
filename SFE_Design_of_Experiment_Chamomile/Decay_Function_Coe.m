function [Upsilon] = Decay_Function_Coe(Re, F, parameters)

    a = parameters{47};
    b = parameters{48};
    c = parameters{49};

    Upsilon  = a + b .* Re + c .* (F * 10^5);
    %Upsilon  = 4.717 + 11.012 .* Re - 0.876 .* (F * 10^5);

end
function [FLUID_SAT] = Saturation_Concentration(SOLID, C_SAT)

    k = 100;            % Growth rate
    %L = 0.5*C_SAT;          % Location of midpoint
    FLUID_SAT = 1 ./ (1 + exp( k .* ( SOLID - C_SAT ) ) );

end


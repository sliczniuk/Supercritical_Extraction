function [FLUID_SAT] = Saturation_Concentration(FLUID, C_SAT)

    k = 4;            % Growth rate
    %L = 0.5*C_SAT;          % Location of midpoint
    FLUID_SAT = 1 ./ (1 + exp( k .* ( FLUID - C_SAT ) ) );

end


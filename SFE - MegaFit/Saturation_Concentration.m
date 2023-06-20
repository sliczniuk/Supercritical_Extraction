function [FLUID_SAT] = Saturation_Concentration(SOLID, C_SAT, km)

    k = 25;            % Growth rate
    %L = 0.5*C_SAT;          % Location of midpoint
    FLUID_SAT = km ./ (1 + exp( k .* ( SOLID - C_SAT ) ) );

end
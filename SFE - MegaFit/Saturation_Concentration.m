function [FLUID_SAT] = Saturation_Concentration(SOLID, shape, km)

    %k = shape;                         % Growth rate
    %L = 0.5*C_SAT;          % Location of midpoint
    %FLUID_SAT = km ./ (1 + exp( shape .* ( SOLID - C_SAT ) ) );
    FLUID_SAT = km .* exp( -shape .* ( SOLID ) );
end
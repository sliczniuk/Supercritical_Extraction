function [Di] = DI_correlation(RHO, k1, k2)

    Di = k1 .* RHO + k2;

end
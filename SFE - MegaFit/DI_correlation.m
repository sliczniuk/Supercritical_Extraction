function [Di] = DI_correlation(RHO, k)

    Di = k(1) .* RHO + k(2);

end
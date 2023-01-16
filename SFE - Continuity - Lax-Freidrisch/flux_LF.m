function [F] = flux_LF(X, Y, dt, dz)

    x_L = X(1:end-2);
    x   = X(2:end-1);
    x_R = X(3:end);

    f_L = Y(1:end-2);
    f   = Y(2:end-1);
    f_R = Y(3:end);

    F_right = Lax_Friedrisch(x  , x_R, f  , f_R, dt, dz);
    F_left  = Lax_Friedrisch(x_L, x  , f_L, f  , dt, dz);
    F       = ( F_right - F_left ) ./ dz;

end
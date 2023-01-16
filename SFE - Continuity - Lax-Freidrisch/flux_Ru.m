function [F] = flux_Ru(X, Y, U, C, epsi, dz)

    x_L = X(1:end-2);
    x   = X(2:end-1);
    x_R = X(3:end);

    f_L = Y(1:end-2);
    f   = Y(2:end-1);
    f_R = Y(3:end);

    u_L = U(1:end-2);
    u   = U(2:end-1);
    u_R = U(3:end);

    q_L = C(1:end-2);
    q   = C(2:end-1);
    q_R = C(3:end);

    e_L = epsi(1:end-2);
    e   = epsi(2:end-1);
    e_R = epsi(3:end);

    F_right = Rusanov(x  , x_R, f  , f_R, u  , u_R, q  , q_R, e_L, e   );
    F_left  = Rusanov(x_L, x  , f_L, f  , u_L, u  , q_L, q  , e  , e_R );
    F       = ( F_right - F_left ) ./ dz;

end
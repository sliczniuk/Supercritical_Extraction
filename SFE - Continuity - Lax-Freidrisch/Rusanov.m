function [F] = Rusanov(x_L, x_R, f_L, f_R, u_L, u_R, q_L, q_R, e_L, e_R)

    lambda = max( q_L + abs(u_L), q_R + abs(u_R) );

    phi    = max( e_L, e_R );

    F = 0.5 .* ( f_L + f_R ) - lambda .* phi .* (x_R - x_L);

end
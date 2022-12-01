function F = RusanovFlux(F_L, F_R, c_L, c_R, u_L, u_R, epsi_L, epsi_R, x_L, x_R)

    v_L = c_L + abs(u_L);
    v_R = c_R + abs(u_R);

    F = (F_L + F_R)./2 - (v_L + v_R)./2 .* (epsi_L + epsi_R)./2 .* (x_R - x_L) ;
end
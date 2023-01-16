function [F] = Lax_Friedrisch(x_L, x_R, f_L, f_R, dt, dz)

    F = (dz/(2*dt)) .* (x_L - x_R) + 0.5 .* ( f_L + f_R );

end
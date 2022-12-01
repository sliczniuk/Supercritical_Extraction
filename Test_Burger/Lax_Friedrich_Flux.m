function F = Lax_Friedrich_Flux(u_i, u_i_1, dx, dt)

    F = 1./2.*( f(u_i_1) + f(u_i) ) - dx./(2.*dt) .* ( u_i - u_i_1 );

end
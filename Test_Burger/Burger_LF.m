function dudt = Burger_LF(t,u,dx,dt)
    
    dudt = zeros(numel(u),1);
    
    F_minus = Lax_Friedrich_Flux(u,[u(2:end); u(end) ],dx,dt);
    
    F_plus  = Lax_Friedrich_Flux([u(1); u(1:end-1) ],u,dx,dt);
    
    dudt = - ( F_minus - F_plus ) ./ dx;

    %dudt = -1/dx .*( ...
%        (1/2) .* (f(u) + f([u(2:end); u(end) ])) - dx/(2*dt) .* ( [u(2:end); u(end) ] - u ) -...
%        (1/2) .* (f([u(1); u(1:end-1) ]) + f(u)) + dx/(2*dt) .* ( u - [u(1); u(1:end-1) ] ) ...
%        );

end
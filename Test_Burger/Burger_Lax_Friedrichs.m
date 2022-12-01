function dudt = Burger_Lax_Friedrichs(t,u)
    dudt = zeros(numel(u),1);
    dudt = -(f(u) - f( [u(1); u(1:end-1)] ) );
    
    unew(2:end-1) = 0.5*(u(3:end)+u(1:end-2)) - 0.5*dt/dx * ...
                    (f(u(3:end)) - f(u(1:end-2)));
    unew(1)   = u(1);
    unew(end) = u(end);
end
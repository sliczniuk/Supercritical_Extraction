function dudt = Burger_E_Conservative(t,u,dx,dt)
    
dudt = zeros(numel(t),numel(u));
    
    for i = 1:numel(t)
    
        unew = u - dt/dx .* (f(u) - f( [u(1); u(1:end-1)] ) );
    
        u = unew;
        dudt(i,:) = u(:);
    
    end

end
function dudt = Burger_E_LF(t,u,dx,dt)
    
dudt = zeros(numel(t),numel(u));
    
    for i = 1:numel(t)
    
        unew =  (1/2)       .* (   [u(2:end); u(end) ]   +    [u(1); u(1:end-1) ] ) - ...
                (1/2)*dt/dx .* (f( [u(2:end); u(end) ] ) - f( [u(1); u(1:end-1) ] ) );
    
        u = unew;
        dudt(i,:) = u(:);
    
    end

end
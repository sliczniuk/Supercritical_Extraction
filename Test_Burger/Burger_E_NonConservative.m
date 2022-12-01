function dudt = Burger_E_NonConservative(t,u,dx,dt)
    
dudt = zeros(numel(t),numel(u));
    
    for i = 1:numel(t)
    
        unew = u - dt/dx .* u .* (u - [u(1); u(1:end-1) ]);
    
        u = unew;
        dudt(i,:) = u(:);
    
    end

end
function xout = SimulateSystemEuler(f,x0,p,dt)

    N = size(p,1);
    Nx = size(x0,1);
    xout = zeros(Nx,N+1);  
    xout(:,1) = x0;

    for k = 1:N
        xdot        = f(xout(:,k), p(k,:)); 
        xout(:,k+1) = full(xout(:,k) + xdot*dt) ;
    end

end
function dudt = burger_Backward(t, u, Nx, dz, a)
    
    dudt = zeros(Nx,1);
    
    dudt(1)    = a * (u(1) - u(Nx)) / dz ;
    dudt(2:Nx) = a .* (u(2:Nx) - u(1:Nx-1)) ./ dz ;

end
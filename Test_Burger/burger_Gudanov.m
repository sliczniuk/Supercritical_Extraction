function dudt = burger_Gudanov(t, u, Nx, dz, a)
    
    dudt = zeros(Nx,1);
    
    dudt(1)    = a * (u(1) - u(Nx)) / dz ;
    dudt(2:Nx-1)    = ( Gudanov_Flux_plus( a,u(2:Nx-1),u(3:Nx) ) - Gudanov_Flux_minus( a,u(2:Nx-1),u(1:Nx-2) ) ) / dz ;
    dudt(Nx) = a .* (u(Nx) - u(Nx-1)) ./ dz ;

end
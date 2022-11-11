function F = Gudanov_Flux_plus(f,u_i,u_plus)

    F = f./2 * ( u_plus + u_i ) - abs(f)./2 .* (u_plus - u_i);

end
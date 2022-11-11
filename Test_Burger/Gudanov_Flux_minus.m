function F = Gudanov_Flux_minus(f,u_i,u_minus)

    F = f./2 * ( u_i + u_minus ) - abs(f)./2 .* (u_i - u_minus);

end
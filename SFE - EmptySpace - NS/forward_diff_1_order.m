function d = forward_diff_1_order(U, U_0, U_B, dz)

    d = [U(2:end); U_B] - U;

    d = d ./ dz;

end
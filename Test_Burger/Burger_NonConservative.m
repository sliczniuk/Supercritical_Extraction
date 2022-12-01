function dudt = Burger_NonConservative(t,u,dx)
    dudt = zeros(numel(u),1);
    dudt = -u .* (u - [u(1); u(1:end-1) ]) ./ dx;
end
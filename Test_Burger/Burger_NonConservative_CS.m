function dudt = Burger_NonConservative_CS(t,u,dx)
    dudt = zeros(numel(u),1);
    dudt = -u .* ( [u(2:end); u(end) ] - [u(1); u(1:end-1) ]) ./ (2*dx);
end
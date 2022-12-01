function dudt = Burger_Conservative(t,u,dx)
    dudt = zeros(numel(u),1);
    dudt = -(f(u) - f( [u(1); u(1:end-1)] ) ) ./ dx;
end
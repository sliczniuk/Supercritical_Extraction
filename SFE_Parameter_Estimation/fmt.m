function output = fmt(X)
    if X == 0 || ~isfinite(X)
       output = '0.000000E+0000';
       return
    end
    D = floor(log10(abs(X)));
    S = X./10.^D;
    output = sprintf('%.6fE%+05d', S, D);
end
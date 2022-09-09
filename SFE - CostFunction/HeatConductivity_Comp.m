function k = HeatConductivity_Comp(T,P,Z,rho, theta)

A1 = theta{28}; % -105.161;
A2 = theta{29}; % 0.9007;
A3 = theta{30}; % 0.0007;
A4 = theta{31}; % 3.50E-15;
A5 = theta{32}; % 3.76E-10;
A6 = theta{33}; % 0.7500;
A7 = theta{34}; % 0.0017;

k = ( A1 + A2.*rho + A3.*(rho.^2) + A4.*(rho.^3).*(T.^3) + A5.*(rho.^4) + A6.*T + A7.*(T.^2) ) ./ (T.^(1/2));

end
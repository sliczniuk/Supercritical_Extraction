function mu = Viscosity(T,P,theta)

% The temperature ranges of the dataare 310–900K and the pressure from 7.5MPa to 101.5MPa. 
% T[K], P[bar]

A1 = theta{35}; % -1.146067E-1;
A2 = theta{36}; %  6.978380E-7;
A3 = theta{37}; %  3.976765E-10;
A4 = theta{38}; %  6.336120E-2;
A5 = theta{39}; % -1.166119E-2;
A6 = theta{40}; %  7.142596E-4;
A7 = theta{41}; %  6.519333E-6;
A8 = theta{42}; % -3.567559E-1;
A9 = theta{43}; %  3.180473E-2;

mu = ( A1 + A2.*P + A3.*P.^2 + A4.*log(T) + A5.*(log(T)).^2 + A6.*(log(T)).^3 ) ./ ( 1 + A7.*P + A8.*log(T) + A9.*(log(T)).^2 );

end
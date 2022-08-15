function mu = Viscosity(T,P,rho,theta, correlation)

switch correlation

    case 'Amooey'
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

        mu = mu*1e3; % [ 10-6 (Pa s) ]

    case 'Fenghour'

        d11       = 0.4071119e-2 ;
        d21       = 0.7198037e-4 ;
        d64       = 0.2411697e-16 ;
        d81       = 0.2971072e-22 ;
        d82       = -0.1627888e-22 ;

        T_red     = T/251.196 ;
        logT_red  = log(T_red) ;

        rho2     = rho*rho ;
        rho6     = rho2*rho2*rho2 ;
        rho8     = rho6*rho2 ;

        a0        = 0.235156 ;
        a1        = -0.491266 ;
        a2        = 5.211155e-2 ;
        a3        = 5.347906e-2 ;
        a4        = -1.537102e-2 ;

        A_mu      = a0 + a1*logT_red + a2*logT_red^2 + a3*logT_red^3  + a4*logT_red^4 ;
        g_mu      = exp(A_mu) ;

        mu_zero   = 1.00697*sqrt(T)/g_mu ;  % /*zero-density viscosity, uPa.s*/

        mu_excess = d11*rho+ d21*rho2 + d64*rho6/(T_red^3) + d81*rho8 + d82*rho8/T_red; %/* in uPa.s*/

        mu         = mu_zero + mu_excess; % /* uPa.s*/

        %mu         = mu*1.e-6 ;           % /* Pa.s*/

    case 'Laesecke'
        

    otherwise
        disp('No such a correlation')

end
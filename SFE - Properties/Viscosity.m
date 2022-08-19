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

        Na     = 6.022e23;
        sigma  = 0.378421e-9;
        M      = 44.0095e-3;
        Tt     = 216.592;
        rho_tl = 1178.53; % density at the triple point

        rho_r  = rho / rho_tl;
        T_star = T   / 200.760;
        Tr     = T   / Tt;

        a	   = [ 1749.354893188350, -369.069300007128, 5423856.34887691, -2.21283852168356, -269503.247933569, 73145.021531826, 5.34368649509278];

        b      = [ -19.572881, 219.73999, -1015.3226, 2471.0125, -3375.1717, 2491.6597, -787.26086,  14.085455, -0.34664158];
        
        t      = [ 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.5, 5.5];

        gamma  = 8.06282737481277;
        c      = [0.360603235428487, 0.121550806591497 ];

        mu_0   = 1.0055*sqrt(T) / ( a(1) + a(2)*T.^(1/6) + a(3)*exp(a(4)*T.^(1/3)) + ( (a(5) + a(6)*T.^(1/3) ) /exp(T.^(1/3)) ) + a(7)*sqrt(T) );

        B      = 0;
        for i=1:8
            B = B + b(i+1) / (T_star^t(i));
        end
        B = b(1) + B;
        
        mu_l   = mu_0 * B * sigma^3 * Na / M; 

        mu_tl  = 0.09436; % rho_tl^(2/3) * sqrt(kB*Tt)/( (M/Na)^(1/6) );

        mu_r   = mu_tl *(  c(1)*Tr*(rho_r^3) + (rho_r^2 + rho_r^gamma)/(Tr - c(2)) );

        mu = mu_0 + rho*mu_l + mu_r;

        mu = mu*1e3;

    otherwise
        disp('No such a correlation')

end